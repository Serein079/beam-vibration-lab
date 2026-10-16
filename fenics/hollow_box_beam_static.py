"""Static deflection of a 3D hollow box beam under gravity and surface load."""

from __future__ import print_function

from dolfin import (
    Constant,
    DirichletBC,
    Expression,
    Function,
    Identity,
    TensorFunctionSpace,
    TestFunction,
    TrialFunction,
    VectorFunctionSpace,
    XDMFFile,
    dot,
    dx,
    ds,
    inner,
    near,
    parameters,
    solve,
    tr,
)
from mshr import Box, Point, generate_mesh
from ufl import nabla_grad


LENGTH = 1.8288
WIDTH = 0.254
HEIGHT = 0.4572

YOUNGS_MODULUS = 200e9
POISSON_RATIO = 0.30
DENSITY = 7_800.0
GRAVITY = -9.8
WALL_THICKNESS = 0.0508
TOP_TRACTION = 581_251.1625


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
    return near(x[0], 0.0) and on_boundary


def build_hollow_beam_mesh():
    outer = Box(Point(0, 0, 0), Point(LENGTH, WIDTH, HEIGHT))
    # One central void yields a closed box section. The upstream two side
    # voids described an I-section and one box had reversed coordinates.
    inner_void = Box(Point(0, WALL_THICKNESS, WALL_THICKNESS),
                     Point(LENGTH, WIDTH - WALL_THICKNESS, HEIGHT - WALL_THICKNESS))
    geometry = outer - inner_void
    return generate_mesh(geometry, 64)


def main():
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["optimize"] = True

    mesh = build_hollow_beam_mesh()
    XDMFFile("hollow_box_beam_mesh.xdmf").write(mesh)

    displacement_space = VectorFunctionSpace(mesh, "P", 1)
    TensorFunctionSpace(mesh, "DG", 0)

    boundary_condition = DirichletBC(
        displacement_space, Constant((0.0, 0.0, 0.0)), clamped_left_boundary
    )

    # dx expects force per volume [N/m^3], not total weight [N].
    gravity_force = DENSITY * GRAVITY
    body_force = Constant((0.0, 0.0, gravity_force))
    traction = Expression(
        (
            "0.0",
            "0.0",
            "near(x[2], height) && x[0] >= 0.97222222 * length ? -traction : 0.0",
        ),
        degree=1,
        height=HEIGHT,
        length=LENGTH,
        traction=TOP_TRACTION,
    )

    trial_displacement = TrialFunction(displacement_space)
    test_displacement = TestFunction(displacement_space)

    stiffness = inner(stress(trial_displacement), strain(test_displacement)) * dx
    load = dot(body_force, test_displacement) * dx + dot(traction, test_displacement) * ds

    displacement = Function(displacement_space, name="Displacement")
    solve(stiffness == load, displacement, boundary_condition)

    output = XDMFFile("hollow_box_beam_static.xdmf")
    output.write(displacement)
    print("Wrote hollow_box_beam_static.xdmf")


if __name__ == "__main__":
    main()
