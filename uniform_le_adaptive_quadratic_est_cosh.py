__author__ = 'lukas'

from triangularmesh.mesh import *
from triangularmesh.utils import *
from dolfin import *
import math
import matplotlib.pyplot as plt

class ExactSolution(Expression):
    def eval(self,values, x):
        values[0] = (math.cosh(10*(1-x[0]))+ math.cosh(10*(1-x[1])))/(2* math.cosh(10))
        return values

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
mesh = TriangularMesh.create_uniform_mesh(fd, [0,1,0,1], .1, pv)

p = mesh.node_list
t = mesh.face_list

maxerr = 0.2

errors_old = []

output = file("convergence/uniform_le_adaptive_quadratic_cosh.dat", 'w')


for i in range(30):
    fmesh = mesh.to_dolfin_mesh()
    file = File('meshes/uniform_le_adaptive_quadratic_cosh_'+str(i)+'.pvd')
    file << fmesh

    V = FunctionSpace(fmesh, 'Lagrange', 1)
    V2 = FunctionSpace(fmesh, 'Lagrange', 2)
    Ve = FunctionSpace(fmesh, 'Lagrange', 5)


    def u0_boundary(x, on_boundary):
        return on_boundary

    u0 = ExactSolution()
    u_e = u0
    bc = DirichletBC(V, u0, u0_boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(0)
    a = inner(grad(u), grad(v))*dx()+100*u*v*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)


    u2 = TrialFunction(V2)
    v = TestFunction(V2)
    a = inner(nabla_grad(u2), nabla_grad(v))*dx
    L = f*v*dx
    u2 = Function(V2)
    solve(a == L, u2, bc)

    l2 = errornorm(u_e, u, norm_type='L2', degree_rise=3)
    ene = errornorm(u_e, u, norm_type='H10', degree_rise=3)
    print "L2 error = ", l2
    print "Energy error = ", ene

    output.write(str(len(mesh.get_nodes()))+" "+str(len(mesh.get_faces()))+" "+str(ene)+"\n")

    if ene < 0.2:
        print len(mesh.get_nodes()), len(mesh.get_faces()), ene

    if ene < 0.05:
        break


    errors_now = compute_energy_errors_on_cells(u2, u, Ve, fmesh)

    if i == 0:
        mesh.refine_all()
        p, t = mesh.get_nodes(), mesh.get_faces()
        #dm.simpplot(p, t, annotate="")
        #plt.show()
    else:
       # print errors_old
        ##print errors_now
        p, t = mesh.get_nodes(), mesh.get_faces()
        #dm.simpplot(p, t, annotate="")
        #plt.show()
        mesh.adaptive_refine(errors_old, errors_now, ratio=0.4)
        #mesh.refine_all_by_longest_edge()


    errors_old = errors_now

print "#nodes = ", len(mesh.get_nodes()), "#triangles = " , len(mesh.get_faces())
output.close()

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
