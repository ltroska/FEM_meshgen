__author__ = 'lukas'

from mesh import *
from dolfin import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


from dolfin import *

# Create mesh and define function space
mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG',2)

P = Expression("0")
Q = Expression("0")

def left_boundary(x, on_boundary):
    return on_boundary and x[0] < DOLFIN_EPS

def bottom_boundary(x, on_boundary):
    return on_boundary and (x[1] < DOLFIN_EPS)

def right_boundary1(x, on_boundary):
    return on_boundary and x[0] > 1 - DOLFIN_EPS and x[1] < 0.5

def right_boundary2(x, on_boundary):
    return on_boundary and x[0] > 1-DOLFIN_EPS and x[1] >= 0.5

def top_boundary1(x, on_boundary):
    return on_boundary and x[1] > 1- DOLFIN_EPS and x[0] < 0.5

def top_boundary2(x, on_boundary):
    return on_boundary and x[1] > 1- DOLFIN_EPS and x[0] > 0.5

def getchi():

    u0_left_boundary = Constant(0)
    u0_bottom_boundary = Expression("x[0]")
    u0_right_boundary1 = Constant(1)
    u0_right_boundary2 =Expression("1.5-x[1]")
    u0_top_boundary1 = Expression("x[0]")
    u0_top_boundary2 = Constant(0.5)

    B1 = DirichletBC(V, u0_left_boundary, left_boundary)
    B2 = DirichletBC(V, u0_bottom_boundary, bottom_boundary)
    B3 = DirichletBC(V, u0_right_boundary1, right_boundary1)
    B4 = DirichletBC(V,u0_right_boundary2, right_boundary2)
    B5 = DirichletBC(V, u0_top_boundary1, top_boundary1)
    B6 = DirichletBC(V,  u0_top_boundary2, top_boundary2)

    BCs = [B1, B2, B3,B4, B5,  B6]

    # Define variational problem
    chi = TrialFunction(V)
    v = TestFunction(V)
    f = P
    a = -1*inner(nabla_grad(chi), nabla_grad(v))*dx
    L = f*v*dx

    # Compute solution
    chi = Function(V)
    solve(a == L, chi, BCs)

    # Plot solution and mesh
   # plot(chi)
    n = V.dim()
    d = mesh.geometry().dim()

    dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]

    return np.array(zip(dof_x, dof_y)), chi.vector().array()
    # Hold plot

def getnu():

    u0_left_boundary = Expression("x[1]")
    u0_bottom_boundary = Constant(0)
    u0_right_boundary1 = Expression("x[1]")
    u0_right_boundary2 = Constant(0.5)
    u0_top_boundary1 = Constant(1)
    u0_top_boundary2 = Expression("1.5-x[0]")

    B1 = DirichletBC(V, u0_left_boundary, left_boundary)
    B2 = DirichletBC(V, u0_bottom_boundary, bottom_boundary)
    B3 = DirichletBC(V, u0_right_boundary1, right_boundary1)
    B4 = DirichletBC(V,u0_right_boundary2, right_boundary2)
    B5 = DirichletBC(V, u0_top_boundary1, top_boundary1)
    B6 = DirichletBC(V,  u0_top_boundary2, top_boundary2)

    BCs = [B1, B2, B3,B4, B5,  B6]

    # Define variational problem
    nu = TrialFunction(V)
    v = TestFunction(V)
    f = Q
    a = inner(nabla_grad(nu), nabla_grad(v))*dx
    L = f*v*dx

    # Compute solution
    nu = Function(V)
    solve(a == L, nu, BCs)

    # Plot solution and mesh

    # Dump solution to file in VTK format
    file = File('poisson.pvd')
    file << nu
   # interactive()

    n = V.dim()
    d = mesh.geometry().dim()

    dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]

    return np.array(zip(dof_x, dof_y)), nu.vector().array()
    # Hold plot

_, x = getchi()
_, y = getnu()

plt.scatter(x, y)
plt.show()