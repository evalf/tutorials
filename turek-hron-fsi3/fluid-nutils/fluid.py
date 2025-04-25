# Turek benchmark fluid solver, modified from https://examples.nutils.org/official-turek/

from precice import Participant
from nutils import cli, export, function
from nutils.mesh import gmsh
from nutils.solver import System
from nutils.SI import Length, Density, Viscosity, Velocity, Time, Stress, Acceleration
from nutils.expression_v2 import Namespace
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
import treelog as log
import numpy


# reference lenghts used in communication
precice_length = Length('m')
precice_time = Time('s')
precice_stress = Stress('Pa')


@dataclass
class Domain:
    channel_length: Length = Length('2.5m')
    channel_height: Length = Length('.41m')
    x_center: Length = Length('.2m')
    y_center: Length = Length('.2m')
    cylinder_radius: Length = Length('5cm')
    structure_length: Length = Length('35cm')
    structure_thickness: Length = Length('2cm')
    min_elemsize: Length = Length('4mm')
    max_elemsize: Length = Length('4cm')

    def generate_mesh(self):
        u = Length('m')  # temporary reference length

        topo, geom = gmsh(Path(__file__).parent / 'fluid.geo', dimension=2, order=2, numbers={
            'channel_length': self.channel_length / u,
            'channel_height': self.channel_height / u,
            'x_center': self.x_center / u,
            'y_center': self.y_center / u,
            'cylinder_radius': self.cylinder_radius / u,
            'structure_length': self.structure_length / u,
            'structure_thickness': self.structure_thickness / u,
            'min_elemsize': self.min_elemsize / u,
            'max_elemsize': self.max_elemsize / u})

        # consistency check
        # zeros = topo.boundary['inlet,wall,outlet,cylinder,structure'].integrate(function.normal(geom) * function.J(geom, 1), degree=2)
        # numpy.testing.assert_allclose(zeros, 0, atol=1e-14, err_msg='boundaries do not form a watertight hull')

        bezier = topo.sample('bezier', 2)
        bezier_structure = topo['fluid'].boundary['structure'].sample(
            'bezier', 3)
        bezier_cylinder = topo['fluid'].boundary['cylinder'].sample(
            'bezier', 3)
        A = topo.points['A'].sample('gauss', 1).eval(geom)
        with export.mplfigure('mesh.jpg', dpi=150) as fig:
            ax = fig.add_subplot(111)
            export.triplot(ax, bezier.eval(geom), hull=bezier.hull)
            export.triplot(ax, bezier_structure.eval(
                geom), hull=bezier_structure.tri, linewidth=1, linecolor='r')
            export.triplot(ax, bezier_cylinder.eval(
                geom), hull=bezier_cylinder.tri, linewidth=1, linecolor='b')
            ax.set_xlim(0, 2 * self.channel_height / u)

        return topo, geom * u


@dataclass
class Fluid:
    density: Density = Density('1kg/L')
    viscosity: Viscosity = Viscosity('1Pa*s')
    velocity: Velocity = Velocity('1m/s')

    def reynolds(self, reference_length):
        return self.density * self.velocity * reference_length / self.viscosity


@dataclass
class Dynamic:
    init: Time = Time('2s')
    window: Time = Time('1s')
    gamma: float = .8
    beta: float = .4

    def set_timestep(self, timestep):
        self.timestep = timestep
        self.timeseries = defaultdict(
            deque(maxlen=round(self.window / timestep)).copy)

    @property
    def times(self):
        'Return all configured time steps for the simulation.'

        t = Time('0s')
        while True:
            t += self.timestep
            yield t

    def ramp_up(self, t):
        return .5 - .5 * numpy.cos(numpy.pi * min(t / self.init, 1))

    def rate(self, v1, subs):
        dt = function.field('dt') * precice_time
        v0 = function.replace_arguments(v1, subs)
        return (v1 - v0) / dt

    def add_and_plot(self, name, t, v, ax):
        'Add data point and plot time series for past window.'

        d = self.timeseries[name]
        d.append((t, v))
        times, values = numpy.stack(d, axis=1)
        ax.plot(times, values)
        ax.set_ylabel(name)
        ax.grid()
        ax.autoscale(enable=True, axis='x', tight=True)

    # The Newmark-beta scheme is used for time integration, taking displacement
    # (solid) or velocity (fluid) as primary variable with time derivatives
    # introduced via helper arguments that are updated after every solve.
    #
    # d = d0 + δt u0 + .5 δt^2 aβ, where aβ = (1-2β) a0 + 2β a
    # => δd = δt u0 + δt^2 [ .5 a0 + β δa ]
    # => δa = [ δd / δt^2 - u0 / δt - .5 a0 ] / β
    #
    # u = u0 + δt aγ, where aγ = (1-γ) a0 + γ a
    # => δu = δt [ a0 + γ δa ]
    # => δa = [ δu / δt - a0 ] / γ

    def newmark_defo_args(self, d, d0=0., u0δt=0., a0δt2=0., **args):
        δaδt2 = (d - d0 - u0δt - .5 * a0δt2) / self.beta
        uδt = u0δt + a0δt2 + self.gamma * δaδt2
        aδt2 = a0δt2 + δaδt2
        return dict(args, d=d + uδt + .5 * aδt2, d0=d, u0δt=uδt, a0δt2=aδt2)

    def newmark_defo(self, d):
        D = self.newmark_defo_args(
            d, *[function.replace_arguments(d, [('d', t)]) for t in ('d0', 'u0δt', 'a0δt2')])
        return D['u0δt'] / self.timestep, D['a0δt2'] / self.timestep**2

    def newmark_velo_args(self, u, u0=0., a0δt=0., **args):
        aδt = a0δt + (u - u0 - a0δt) / self.gamma
        return dict(args, u=u + aδt, u0=u, a0δt=aδt)

    def newmark_velo(self, u):
        D = self.newmark_velo_args(
            u, *[function.replace_arguments(u, [('u', t)]) for t in ('u0', 'a0δt')])
        return D['a0δt'] / self.timestep


def main(domain: Domain = Domain(), fluid: Fluid = Fluid(), dynamic: Dynamic = Dynamic()):

    participant = Participant('Fluid', '../precice-config.xml', 0, 1)

    log.info('Re:', fluid.reynolds(2 * domain.cylinder_radius))

    topo, geom = domain.generate_mesh()

    d = topo.field('d', btype='std', degree=1, shape=[2]) * precice_length

    sqr = topo.boundary['structure'].sample(
        'bezier', 2).integral((d - geom) @ (d - geom)) / 'm2'
    r_points = System(sqr, trial='d').solve_constraints(droptol=1e-10)['d']
    r_where = numpy.isfinite(r_points).any(1)
    assert not numpy.isnan(r_points[r_where]).any()
    r_name = "Fluid-Mesh-Nodes"
    r_ids = participant.set_mesh_vertices(r_name, r_points[r_where])

    w_name = "Fluid-Mesh-Centers"
    w_sample = topo.boundary['structure'].sample('gauss', degree=1)
    w_ids = participant.set_mesh_vertices(
        w_name, w_sample.eval(geom) / precice_length)

    participant.initialize()

    timestep = participant.get_max_time_step_size()
    dynamic.set_timestep(timestep * precice_time)

    ns = Namespace()
    ns.δ = function.eye(2)
    ns.ρf = fluid.density
    ns.μf = fluid.viscosity

    ns.xref = geom
    ns.define_for('xref', gradient='∇ref', jacobians=('dVref', 'dSref'))

    ns.d = d
    ns.x_i = 'xref_i + d_i'
    ns.F_ij = '∇ref_j(x_i)'  # deformation gradient tensor
    ns.C_ij = 'F_ki F_kj'  # right Cauchy-Green deformation tensor
    ns.J = numpy.linalg.det(ns.F)

    ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))

    ns.urel = topo.field('u', btype='std', degree=2,
                         shape=(2,)) * fluid.velocity
    ns.v, ns.a = dynamic.newmark_defo(ns.d)
    ns.arel = dynamic.newmark_velo(ns.urel)
    ns.u_i = 'v_i + urel_i'
    ns.DuDt_i = 'a_i + arel_i + ∇_j(u_i) urel_j'  # material derivative

    ns.p = topo.field('p', btype='std', degree=1) * \
        fluid.viscosity * fluid.velocity / domain.cylinder_radius
    ns.σ_ij = 'μf (∇_j(u_i) + ∇_i(u_j)) - p δ_ij'  # fluid stress tensor

    y = ns.xref[1] / domain.channel_height
    uin = 6 * fluid.velocity * y * (1 - y)
    sqr = topo.boundary['wall,cylinder,structure'].integral(
        ns.urel @ ns.urel, degree=4) / 'm2/s2'
    sqr += topo.boundary['wall,cylinder,inlet,outlet'].integral(
        ns.d @ ns.d, degree=4) / 'm2'
    sqr += topo.boundary['inlet'].integral(
        (ns.urel[0] - uin)**2 * ns.dSref, degree=4) / 'm3/s2'
    sqr += topo.boundary['inlet,outlet'].integral(
        ns.urel[1]**2, degree=4) / 'm2/s2'
    cons = System(sqr, trial='u,d').solve_constraints(droptol=1e-10)
    ucons = cons['u']
    cons['u'] = ucons * 0

    # ρf DuDt = div σ =>  ∀q: 0 = ∫ q·(ρf DuDt - div σ) = ∫ (ρf q·DuDt + ∇q:σ)
    ns.utest = function.replace_arguments(
        ns.urel, 'u:utest') / (fluid.viscosity * fluid.velocity**2)
    ns.ptest = function.replace_arguments(
        ns.p, 'p:ptest') / (fluid.viscosity * fluid.velocity**2)
    res = topo.integral(
        '(utest_i ρf DuDt_i + ∇_j(utest_i) σ_ij + ptest ∇_k(u_k)) dV' @ ns, degree=4)
    sys_up = System(res, trial='u,p', test='utest,ptest')

    # t = -σ·n => ∀q: ∮ -q·t = ∮ q·σ·n = ∫ div(q·σ) = ∫ (∇q:σ + q·div σ) = ∫ (∇q:σ + ρf q·DuDt)
    ns.t = topo.field('t', btype='std', degree=2, shape=[
                      2]) * (fluid.viscosity * fluid.velocity / domain.cylinder_radius)
    res += topo.boundary['cylinder,structure'].integral(
        'utest_i t_i dS' @ ns, degree=4)
    sys_t = System(res, trial='t', test='utest')
    sqr = topo.boundary['cylinder,structure'].integral(
        't_i t_i' @ ns, degree=4) / 'Pa2'
    cons['t'] = numpy.isnan(
        System(sqr, trial='t').solve_constraints(droptol=1e-10)['t'])
    F = topo.boundary['cylinder,structure'].integral('t_i dS' @ ns, degree=4)

    # mesh continuation with jacobian based stiffness
    sqr = topo.integral('C_kk - 2 log(J)' @ ns, degree=4)
    sys_d = System(sqr, trial='d')

    # initial values
    args = {a: numpy.zeros(function.arguments_for(res)[a].shape) for a in 'ud'}

    bezier = topo['fluid'].sample('bezier', 3)
    bezier = bezier.subset(bezier.eval(geom[0]) < 2.2 * domain.channel_height)
    x_bz = bezier.bind(ns.x)
    u_bz = bezier.bind(ns.u)
    p_bz = bezier.bind(ns.p) - \
        topo.points['B'].sample('gauss', 1).bind(ns.p)[0]

    bbezier = topo['fluid'].boundary['cylinder,structure'].sample('bezier', 3)
    x_bbz = bbezier.bind(ns.x)

    assert participant.requires_writing_checkpoint()
    checkpoint = args

    with log.iter.plain('timestep', dynamic.times) as times:
        t = next(times)

        while participant.is_coupling_ongoing():

            if participant.requires_reading_checkpoint():
                args = checkpoint

            cons['d'][r_where] = participant.read_data(
                r_name, 'Displacement', r_ids, timestep)
            cons['u'] = ucons * dynamic.ramp_up(t)

            args = dynamic.newmark_defo_args(**args)
            args = dynamic.newmark_velo_args(**args)
            args = sys_d.solve(arguments=args, constrain=cons,
                               tol=1e-10)  # mesh extension
            # Navier-Stokes time step
            args = sys_up.solve(arguments=args, constrain=cons, tol=1e-10)
            args = sys_t.solve(arguments=args, constrain=cons,
                               tol=1e-10)  # project traction

            participant.write_data(
                w_name, 'Stress', w_ids, w_sample.eval(ns.t, args) / precice_stress)
            participant.advance(timestep)

            if participant.requires_writing_checkpoint():
                checkpoint = args

                x, xb, u, p = function.eval([x_bz, x_bbz, u_bz, p_bz], args)
                with export.mplfigure('solution.jpg', dpi=150) as fig:
                    pstep = 25 * fluid.viscosity * fluid.velocity / domain.channel_height
                    ax = fig.add_subplot(
                        111, title=f'flow at t={t:.3s}, pressure contours every {pstep:.0Pa}', ylabel='[m]')
                    vmax = 2 * fluid.velocity * dynamic.ramp_up(t)
                    im = export.triplot(ax, x / 'm', numpy.linalg.norm(u / 'm/s', axis=1),
                                        tri=bezier.tri, cmap='inferno', clim=(0, vmax / 'm/s'))
                    ax.tricontour(*(x / 'm').T, bezier.tri, p / pstep, numpy.arange(*numpy.quantile(numpy.ceil(p / pstep), [0, 1])),
                                  colors='white', linestyles='solid', linewidths=1, alpha=.33)
                    fig.colorbar(im, orientation='horizontal',
                                 label=f'velocity [m/s]')
                    export.triplot(ax, xb / 'm', hull=bbezier.tri, linewidth=1)
                    ax.set_xlim(0, 2 * domain.channel_height / 'm')
                    ax.set_ylim(0, domain.channel_height / 'm')

                D, L = function.eval(F, arguments=args)
                log.info(f'lift: {L:N/m}')
                log.info(f'drag: {D:N/m}')
                with export.mplfigure('force.jpg', dpi=150) as fig:
                    dynamic.add_and_plot(
                        'lift [N/m]', t / 's', L / 'N/m', ax=fig.add_subplot(211))
                    dynamic.add_and_plot(
                        'drag [N/m]', t / 's', D / 'N/m', ax=fig.add_subplot(212, xlabel='time [s]'))

                t = next(times)


if __name__ == '__main__':
    cli.run(main)
