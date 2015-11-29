__author__ = 'lukas'

from triangularmesh.mesh import *
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

#mesh = TriangularMesh.create_spring_based_mesh(fd2, [0,1,0,1], .3, element_size_function=dm.huniform)
mesh = TriangularMesh.create_uniform_mesh(fd, [0,1,0,1], .1, pv)

p = mesh.node_list
t = mesh.face_list

maxerr = 0.2
output = file("convergence/uniform_del.dat", 'w')


errors_old = []

for i in range(10):
    mesh.refine_all(refine_method='del')
    p, t = mesh.get_nodes(), mesh.get_faces()

    fmesh = mesh.to_dolfin_mesh()
    file = File('meshes/uniform_del_'+str(i)+'.pvd')
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
    output.write(str(len(mesh.get_nodes()))+" "+str(len(mesh.get_faces()))+" "+str(ene)+"\n")

    if ene < 0.2:
        break

    errors_now, _ = computeErrorsOnCells(u_e, u, Ve, fmesh)


    errors_old = errors_now

print "#nodes = ", len(mesh.get_nodes()), "#triangles = " , len(mesh.get_faces())
print len(t)
# Plot solution and mesh
plot(u)
#plot(interpolate(u_e, Ve))
plot(fmesh)
output.close()

#print errors_now

# Dump solution to file in VTK format
file = File('solutions/uniform_del_poisson.pvd')
file << u

# Hold plot
interactive()
