__author__ = 'lukas'

from utils import *
from mesh import *
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

#mesh = TriangularMesh.create_spring_based_mesh(fd2, [0,1,0,1], .3, element_size_function=dm.huniform)
mesh = TriangularMesh.create_equilateral_uniform_mesh(fd, [0,1,0,1], .1, pv)

p = mesh.node_list
t = mesh.face_list

maxerr = 0.2

errors_old = []

for i in range(30):
    fmesh = mesh.to_dolfin_mesh()
    file = File('meshes/uniform_le_adaptive_'+str(i)+'.pvd')
    file << fmesh

    V = FunctionSpace(fmesh, 'Lagrange', 1)
    Ve = FunctionSpace(fmesh, 'Lagrange', 3)

    # Define boundary conditions
    u0 = Expression('exp(-100*((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))')
    u_e = Expression('exp(-100*((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))')

    def u0_boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, u0, u0_boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression('-40000*exp(-100 * ((x[0]-0.25)*(x[0]-0.25)+(x[1]-0.25)*(x[1]-0.25)))*(x[0]*x[0]-0.5*x[0]+x[1]*x[1]-0.5*x[1]+0.115)')
    a = inner(nabla_grad(u), nabla_grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)
    l2 = errornorm(u_e, u, norm_type='L2', degree_rise=3)
    ene = errornorm(u_e, u, norm_type='H10', degree_rise=3)
    print "L2 error = ", l2
    print "Energy error = ", ene

    plot(fmesh)
    interactive()
    if l2 < 0.001:
        break

    errors_now, _ = computeErrorsOnCells(u_e, u, Ve, fmesh)
    if i == 0:
        mesh.refine_all_with_bigger_error(errors_now, 1e-7, refine_method='del')
        p, t = mesh.get_nodes(), mesh.get_faces()
        #dm.simpplot(p, t, annotate="")
        #plt.show()
    else:
       # print errors_old
        ##print errors_now
        p, t = mesh.get_nodes(), mesh.get_faces()
        #dm.simpplot(p, t, annotate="")
        #plt.show()
        mesh.refine_all_with_bigger_error(errors_now, 1e-4, refine_method='del')

       #mesh.refine_all_by_longest_edge()

    errors_old = errors_now

p, t = mesh.get_nodes(), mesh.get_faces()
print "#nodes = ", len(p), " #triangles = ", len(t)

# Plot solution and mesh
plot(u)
#plot(interpolate(u_e, Ve))
plot(fmesh)

#print errors_now

# Dump solution to file in VTK format
file = File('solutions/uniform_le_adaptive_poisson.pvd')
file << u

# Hold plot
interactive()
