from triangularmesh.mesh import *
from triangularmesh.utils import *
from dolfin import *
import math
import matplotlib.pyplot as plt

mesh = TriangularMesh()

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

fd = lambda p: dm.dpoly(p, pv)
fd2 = lambda p: dm.drectangle0(p, 0,1, 0,1)

mesh = TriangularMesh.create_uniform_mesh(fd, [0,1,0,1], .1, pv).to_dolfin_mesh()
V = FunctionSpace(mesh, "Lagrange", 1)


# Define boundary conditions
u0 = Expression('exp(-100*((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))')
u_e = Expression('exp(-100*((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))')


def u0_boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

u = TrialFunction(V)
v = TestFunction(V)
f = Expression('-40000*exp(-100 * ((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))*(x[0]*x[0]-0.5*x[0]+x[1]*x[1]-0.5*x[1]+0.115)')
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

# Define function for the solution
u = Function(V)

# Define goal functional (quantity of interest)
M = inner(u, u)*dx()

# Define error tolerance
#tol = 0.00005
tol = 0.000082

# Solve equation a = L with respect to u and the given boundary
# conditions, such that the estimated error (measured in M) is less
# than tol
problem = LinearVariationalProblem(a, L, u, bc)
solver = AdaptiveLinearVariationalSolver(problem, M)
solver.parameters["error_control"]["dual_variational_solver"]["linear_solver"] = "cg"
solver.solve(tol)

solver.summary()

print errornorm(u_e, u.leaf_node(), norm_type='H10', degree_rise=3), len(mesh.leaf_node().coordinates())

file = File('meshes/fenics_refine.pvd')
file << mesh.leaf_node()

