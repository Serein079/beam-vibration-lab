"""Static 3D cantilever beam deflection with a distributed surface traction."""

from __future__ import print_function

from fenics import (
    BoxMesh,
    Constant,
    DirichletBC,
    Expression,
    Function,
    Identity,
    Point,
    TensorFunctionSpace,
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
SURFACE_TRACTION = 20_000.0


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


def main():
    mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(LENGTH, WIDTH, HEIGHT), 20, 6, 6)
    displacement_space = VectorFunctionSpace(mesh, "P", 1)
    TensorFunctionSpace(mesh, "P", 1)

    boundary_condition = DirichletBC(
        displacement_space, Constant((0.0, 0.0, 0.0)), clamped_left_boundary
    )

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

    displacement = Function(displacement_space, name="Displacement")
    solve(stiffness == load, displacement, boundary_condition)

    output = XDMFFile("static_cantilever_deflection.xdmf")
    output.write(displacement)
    print("Wrote static_cantilever_deflection.xdmf")


if __name__ == "__main__":
    main()
