from triangularmesh.mesh import *
from triangularmesh.meshquality import *
import matplotlib.pyplot as plt

from math import *

def solve_problem(mesh):
    mesh = mesh.to_dolfin_mesh()


    V = FunctionSpace(mesh, 'Lagrange', 1)

    # Define boundary conditions
    u0 = Expression("0.3*sin(10*(x[0]-0.1)*(x[1]+0.2))")
    u_e = u0

    def u0_boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, u0, u0_boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression("-0.3*100*(-0.1+x[0])*(-0.1+x[0])*sin(10*(-0.1+x[0])*(0.2+x[1]))-0.3*100*(0.2+x[1])*(0.2+x[1])*sin(10*(-0.1+x[0])*(0.2+x[1]))")

    a = -inner(nabla_grad(v), nabla_grad(u))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    l2 = errornorm(u_e, u, norm_type='L2', degree_rise=3)
    ene = errornorm(u_e, u, norm_type='H10', degree_rise=3)
    print "L2 error = ", l2, " Energieerror = ", ene



dp = 0.1
pv = np.array([[0,0]])

for i in range(1,int(math.ceil(1/dp))):
    pv = np.append(pv, [[i*dp, 0]], axis=0)
for i in range(int(math.ceil(.5/dp))):
    pv = np.append(pv, [[1, i*dp]], axis=0)
for i in range(int(math.ceil(.5/dp))):
    pv = np.append(pv, [[1-i*dp, .5]], axis=0)
for i in range(int(math.ceil(.5/dp))):
    pv = np.append(pv, [[.5, .5+i*dp]], axis=0)
for i in range(int(math.ceil(.5/dp))):
    pv = np.append(pv, [[.5-i*dp, 1]], axis=0)
for i in range(int(math.ceil(1/dp))):
    pv = np.append(pv, [[0, 1-i*dp]], axis=0)

p2 = np.array([[1,0.5], [0.5,1], [0.5,0.5]])

fd = lambda p: dm.dpoly(p, pv)

fd2 = lambda p: np.maximum(dm.drectangle0(p, 0,1, 0,1), -dm.drectangle0(p, 0.5,1, 0.5,1))

meshes = []

meshes.append([TriangularMesh.create_equilateral_uniform_mesh(fd, [0, 1, 0, 1],0.075, pv), "uniform"])
meshes.append([TriangularMesh.create_elliptic_mesh(15, Expression("4*x[0]+x[1]"), Expression("x[0]+4*x[1]"), fd), "elliptic"])
meshes.append([TriangularMesh.create_quasi_random_mesh(fd, [0,1,0,1], 185, 2, 3, pv), "quasi-random"])
meshes.append([TriangularMesh.create_spring_based_mesh(fd2, [0,1,0,1], 0.07, pv, dm.huniform), "uniform distmesh"])

fh = lambda p: 7.5*(p[:,0]-0.5)**2+7.5*(p[:,1]-0.5)**2+1
meshes.append([TriangularMesh.create_spring_based_mesh(fd2, [0,1,0,1], 0.035, pv, fh), "corner distmesh"])

for item in meshes:
    mesh, name = item

    if name == "uniform":
        dm.simpplot(mesh.get_nodes(), mesh.get_faces())
        plt.show()
    print name, len(mesh.get_nodes())
    print "minimum angle " , compute_minimum_angle(mesh)
    print "aspect ratio " , compute_aspect_ratio(mesh)
    print "edge ratio " , compute_edge_ratio(mesh)
    print "skewness " , compute_skewness(mesh)
    solve_problem(mesh)
    print "------------------------------------------------------------"

