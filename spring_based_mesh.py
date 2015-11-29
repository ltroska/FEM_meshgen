__author__ = 'lukas'

from triangularmesh.mesh import *
from dolfin import *

mesh = TriangularMesh()

dp = 0.1

pv = np.array([[0,0], [1,0], [1, 0.5], [0.5, 0.5], [0.5, 1], [0, 1]])

fd = lambda p: np.maximum(dm.drectangle0(p, 0,1, 0,1), -dm.drectangle0(p, 0.5,1, 0.5,1))

fh = lambda p: 1-8*fd(p)

fh2 = lambda p: 7.5*(p[:,0]-0.5)**2+7.5*(p[:,1]-0.5)**2+1

maxerr = 0.2
dxr = .025

for i in range(30):
    mesh = TriangularMesh.create_spring_based_mesh(fd, [0,1,0,1], dxr, pv, dm.huniform)
    fmesh = mesh.to_dolfin_mesh()
    file = File('meshes/spring_based_'+str(i)+'.pvd')
    file << fmesh
    file = File('meshes/spring_based_'+str(i)+'.xml')
    file << fmesh
    plot(fmesh)
    interactive()

    V = FunctionSpace(fmesh, 'Lagrange', 1)

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

    if l2 < 0.003:
        break

    dxr /= 2

    #errors_now, _ = computeErrorsOnCells(u_e, u, Ve, fmesh)

print len(mesh.get_faces())
# Plot solution and mesh
plot(u)
#plot(interpolate(u_e, Ve))
plot(fmesh)

#print errors_now

# Dump solution to file in VTK format
#file = File('solutions/spring_based_poisson.pvd')
#file << u

# Hold plot
interactive()
