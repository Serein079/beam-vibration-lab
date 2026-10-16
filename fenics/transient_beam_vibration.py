"""Transient vibration of a cantilever beam after static preloading."""

from __future__ import print_function

from fenics import (
    BoxMesh,
    Constant,
    DirichletBC,
    Expression,
    Function,
    Identity,
    Point,
    TestFunction,
    TrialFunction,
    VectorFunctionSpace,
    XDMFFile,
    dot,
    dx,
    ds,
    inner,
    lhs,
    near,
    rhs,
    solve,
    tr,
)
from ufl import nabla_grad


LENGTH = 50e-3
WIDTH = 5e-3
HEIGHT = 5e-3

YOUNGS_MODULUS = 200e9
POISSON_RATIO = 0.30
DENSITY = 7_800.0
SURFACE_TRACTION = 20_000.0

TIME_FINAL = 2.75e-3
TIME_STEP = 5e-7


def strain(displacement):
    return 0.5 * (nabla_grad(displacement) + nabla_grad(displacement).T)


def stress(displacement):
    lame_lambda = (
        YOUNGS_MODULUS
        * POISSON_RATIO
        / ((1.0 + POISSON_RATIO) * (1.0 - 2.0 * POISSON_RATIO))
    )
    lame_mu = YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO))
    epsilon = strain(displacement)
    return lame_lambda * tr(epsilon) * Identity(3) + lame_mu * (
        epsilon + epsilon.T
    )


def clamped_left_boundary(x, on_boundary):
    return on_boundary and near(x[0], 0.0)


def solve_static_initial_condition(mesh, displacement_space, boundary_condition):
    body_force = Constant((0.0, 0.0, 0.0))
    traction = Expression(
        (
            "0.0",
            "near(x[1], width) && x[0] >= 0.9 * length ? -traction : 0.0",
            "0.0",
        ),
        degree=1,
        length=LENGTH,
        width=WIDTH,
        traction=SURFACE_TRACTION,
    )

    trial_displacement = TrialFunction(displacement_space)
    test_displacement = TestFunction(displacement_space)
    stiffness = inner(stress(trial_displacement), strain(test_displacement)) * dx
    load = dot(body_force, test_displacement) * dx + dot(traction, test_displacement) * ds

    displacement = Function(displacement_space)
    solve(stiffness == load, displacement, boundary_condition)
    return displacement


def main():
    mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(LENGTH, WIDTH, HEIGHT), 20, 6, 6)
    displacement_space = VectorFunctionSpace(mesh, "P", 1)
    boundary_condition = DirichletBC(
        displacement_space, Constant((0.0, 0.0, 0.0)), clamped_left_boundary
    )

    previous_displacement = solve_static_initial_condition(
        mesh, displacement_space, boundary_condition
    )
    two_steps_back = Function(displacement_space)
    two_steps_back.assign(previous_displacement)

    # Taylor initialization for zero initial velocity: u(-dt)=u(0)+dt^2*a(0)/2.
    # a(0) is the acceleration immediately after releasing the static load.
    acceleration_trial = TrialFunction(displacement_space)
    acceleration_test = TestFunction(displacement_space)
    initial_acceleration = Function(displacement_space)
    solve(DENSITY * dot(acceleration_trial, acceleration_test) * dx ==
          -inner(stress(previous_displacement), strain(acceleration_test)) * dx,
          initial_acceleration, boundary_condition)
    two_steps_back.vector().axpy(0.5 * TIME_STEP**2, initial_acceleration.vector())
    two_steps_back.vector().apply("insert")

    current_trial = TrialFunction(displacement_space)
    test_displacement = TestFunction(displacement_space)
    body_force = Constant((0.0, 0.0, 0.0))
    released_traction = Constant((0.0, 0.0, 0.0))

    residual = (
        -(TIME_STEP * TIME_STEP) * inner(stress(current_trial), strain(test_displacement)) * dx
        + (TIME_STEP * TIME_STEP) * dot(body_force, test_displacement) * dx
        - DENSITY * dot(current_trial, test_displacement) * dx
        + 2.0 * DENSITY * dot(previous_displacement, test_displacement) * dx
        - DENSITY * dot(two_steps_back, test_displacement) * dx
        + (TIME_STEP * TIME_STEP) * dot(released_traction, test_displacement) * ds
    )

    system_matrix, system_rhs = lhs(residual), rhs(residual)
    current_displacement = Function(displacement_space, name="Displacement")

    output = XDMFFile("transient_beam_vibration.xdmf")
    output.parameters["rewrite_function_mesh"] = False
    output.parameters["flush_output"] = True

    output.write(previous_displacement, 0.0)
    n_steps = int(round(TIME_FINAL / TIME_STEP))
    for step in range(1, n_steps + 1):
        time = step * TIME_STEP
        print("t = %.8f" % time)
        solve(system_matrix == system_rhs, current_displacement, boundary_condition)
        output.write(current_displacement, time)

        two_steps_back.assign(previous_displacement)
        previous_displacement.assign(current_displacement)

    print("Wrote transient_beam_vibration.xdmf")


if __name__ == "__main__":
    main()
