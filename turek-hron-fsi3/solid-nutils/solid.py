# Turek benchmark solid solver, modified from https://examples.nutils.org/official-turek/

from precice import Participant
from nutils import cli, export, function
from nutils.mesh import gmsh
from nutils.solver import System
from nutils.SI import Length, Density, Viscosity, Velocity, Time, Stress, Acceleration, Pressure
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
    x_center: Length = Length('.2m')
    y_center: Length = Length('.2m')
    cylinder_radius: Length = Length('5cm')
    structure_length: Length = Length('35cm')
    structure_thickness: Length = Length('2cm')
    elemsize: Length = Length('4mm')

    def generate_mesh(self):
        u = Length('m')  # temporary reference length

        topo, geom = gmsh(Path(__file__).parent / 'solid.geo', dimension=2, order=2, numbers={
            'x_center': self.x_center / u,
            'y_center': self.y_center / u,
            'cylinder_radius': self.cylinder_radius / u,
            'structure_length': self.structure_length / u,
            'structure_thickness': self.structure_thickness / u,
            'elemsize': self.elemsize / u})

        bezier = topo.sample('bezier', 2)
        bezier_fluid = topo.boundary['fluid'].sample('bezier', 3)
        bezier_cylinder = topo.boundary['cylinder'].sample('bezier', 3)
        with export.mplfigure('mesh.jpg', dpi=150) as fig:
            ax = fig.add_subplot(111)
            export.triplot(ax, bezier.eval(geom), hull=bezier.hull)
            export.triplot(ax, bezier_fluid.eval(geom),
                           hull=bezier_fluid.tri, linewidth=1, linecolor='r')
            export.triplot(ax, bezier_cylinder.eval(
                geom), hull=bezier_cylinder.tri, linewidth=1, linecolor='b')
            ax.set_xlim(0, 4 * self.x_center / u)
            ax.set_ylim(0, 2 * self.y_center / u)

        return topo, geom * u


@dataclass
class Solid:
    density: Density = Density('10kg/L')
    poisson_ratio: float = .4
    shear_modulus: Pressure = Pressure('.5MPa')
    gravity: Acceleration = Acceleration('0m/s2')

    def lame_parameters(self):
        return 2 * self.shear_modulus * self.poisson_ratio / (1 - 2 * self.poisson_ratio), self.shear_modulus

    def young(self):
        return 2 * self.shear_modulus * (1 + self.poisson_ratio)


@dataclass
class Dynamic:
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


@dataclass
class DummyParticipant:

    timestep: Time = Time('10ms')
    endtime: Time = Time('10s')

    def __post_init__(self):
        self.countdown = self.endtime / self.timestep
        self.npoints = {}

    def initialize(self):
        pass

    def set_mesh_vertices(self, mesh_name, points):
        self.npoints[mesh_name] = len(points)

    def is_coupling_ongoing(self):
        return self.countdown > 0

    def get_max_time_step_size(self):
        return self.timestep / precice_time

    def requires_writing_checkpoint(self):
        return True

    def requires_reading_checkpoint(self):
        return False

    def read_data(self, mesh_name, *_):
        return numpy.zeros((self.npoints[mesh_name], 2))

    def write_data(self, *_):
        pass

    def advance(self, dt):
        self.countdown -= 1


def main(domain: Domain = Domain(), solid: Solid = Solid(), dynamic: Dynamic = Dynamic()):

    participant = Participant('Solid', '../precice-config.xml', 0, 1)

    topo, geom = domain.generate_mesh()

    rw_name = "Solid-Mesh"
    rw_sample = topo.boundary.sample('gauss', degree=1)
    rw_ids = participant.set_mesh_vertices(
        rw_name, rw_sample.eval(geom) / precice_length)

    participant.initialize()

    timestep = participant.get_max_time_step_size()
    dynamic.set_timestep(timestep * precice_time)

    ns = Namespace()
    ns.δ = function.eye(2)
    ns.xref = geom
    ns.define_for('xref', gradient='∇ref', jacobians=('dVref', 'dSref'))

    ns.ρs = solid.density
    ns.λs, ns.μs = solid.lame_parameters()
    ns.g = -solid.gravity * ns.δ[1]

    ns.d = topo.field('d', btype='std', degree=2, shape=(
        2,)) * domain.cylinder_radius  # deformation at the end of the timestep
    ns.v, ns.a = dynamic.newmark_defo(ns.d)

    ns.dtest = function.replace_arguments(
        ns.d, 'd:dtest') / (solid.shear_modulus * domain.cylinder_radius**2)
    ns.traction = function.field(
        'traction', rw_sample.basis(), shape=(2,)) * precice_stress

    ns.x = ns.xref + ns.d
    ns.define_for('x', gradient='∇', jacobians=('dV', 'dS'))

    ns.F_ij = '∇ref_j(x_i)'  # deformation gradient tensor
    ns.C_ij = 'F_ki F_kj'  # right Cauchy-Green deformation tensor
    ns.E_ij = '.5 (C_ij - δ_ij)'  # Green-Lagrangian strain tensor
    ns.S_ij = 'λs E_kk δ_ij + 2 μs E_ij'  # 2nd Piola–Kirchhoff stress
    ns.P_ij = 'F_ik S_kj'  # 1st Piola–Kirchhoff stress

    res = topo.integral(
        '(∇ref_j(dtest_i) P_ij + dtest_i ρs (a_i - g_i)) dVref' @ ns, degree=4)
    res -= rw_sample.integral('dtest_i traction_i dS' @ ns)

    sqr = topo.boundary['cylinder'].integral('d_k d_k' @ ns, degree=4) / 'm2'
    cons = System(sqr, trial='d').solve_constraints(droptol=1e-10)

    system = System(res, trial='d', test='dtest')

    # initial values
    args = {'d': numpy.zeros(function.arguments_for(res)['d'].shape)}

    bbezier = topo.boundary.sample('bezier', 3)
    x_bbz = bbezier.bind(ns.x)

    assert participant.requires_writing_checkpoint()
    checkpoint = args

    with log.iter.plain('timestep', dynamic.times) as times:
        t = next(times)

        while participant.is_coupling_ongoing():

            if participant.requires_reading_checkpoint():
                args = checkpoint

            args = dynamic.newmark_defo_args(**args)
            args['traction'] = participant.read_data(
                rw_name, 'Stress', rw_ids, timestep)
            args = system.solve(arguments=args, constrain=cons, tol=1e-10)

            participant.write_data(
                rw_name, 'Displacement', rw_ids, rw_sample.eval(ns.d, args) / precice_length)
            participant.advance(timestep)

            if participant.requires_writing_checkpoint():
                checkpoint = args

                xb = function.eval(x_bbz, args)
                with export.mplfigure('solution.jpg', dpi=150) as fig:
                    ax = fig.add_subplot(
                        111, title=f'deformation at t={t:.3s}', ylabel='[m]')
                    export.triplot(ax, xb / 'm', hull=bbezier.tri, linewidth=1)
                    ax.set_xlim(0, 4 * domain.x_center / 'm')
                    ax.set_ylim(0, 2 * domain.y_center / 'm')

                ux, uy = topo.points['A'].sample(
                    'gauss', 1).eval(ns.d, arguments=args)[0]
                log.info(f'uy: {uy:mm}')
                log.info(f'ux: {ux:mm}')
                with export.mplfigure('tip-displacement.jpg', dpi=150) as fig:
                    dynamic.add_and_plot(
                        'uy [mm]', t / 's', uy / 'mm', ax=fig.add_subplot(211))
                    dynamic.add_and_plot(
                        'ux [mm]', t / 's', ux / 'mm', ax=fig.add_subplot(212, xlabel='time [s]'))

                t = next(times)


if __name__ == '__main__':
    cli.run(main)
