from triangularmesh.mesh import *
from triangularmesh.meshquality import *
import matplotlib.pyplot as plt

from math import *



def solve_problem(mesh_name, u_e = None, p = 1):
    mesh = Mesh(mesh_name+".xml")
    facets = MeshFunction("size_t", mesh, mesh_name+"_facets.xml")

    parameters['allow_extrapolation'] = True

   # plot(mesh)
   # plot(facets)
   # interactive()
    # An elasticity problem
    VV = VectorFunctionSpace(mesh, 'CG', p)
    trial = TrialFunction(VV)
    test = TestFunction(VV)
    # boundary conditions
    right = CompiledSubDomain('near(x[0],0.) && on_boundary')
    facet = FacetFunction('size_t', mesh, 0)
    right.mark(facet, 1)
    bcs = [DirichletBC(VV, Constant((0,0)), facets, 4), DirichletBC(VV, Constant((0,0)), facets, 6)]
    dd = ds[facet]
    # constants
    poisson = 0.3
    young = 1.
    ll = (young*poisson)/((1.+poisson)*(1.-2.*poisson))
    mu = young/(2*(1.+poisson))
    # forms
    aa = (mu*inner(grad(trial), grad(test))+(mu+ll)*inner(div(trial), div(test)))*dx
    LL = inner(Constant((0,0)),test)*dx+inner(Constant((0.1,0)),test)*dd(1)
    # solve
    uu = Function(VV)
    solve(aa == LL, uu, bcs)

    #plot(uu, mode='displacement', hardcopy_prefix='elasticity', title='elasticity')
   # interactive()

    u_Ve = uu
   # plot(uu, mode='displacement', hardcopy_prefix='elasticity', title='elasticity')
    #interactive()

    error = (mu*inner(grad(u_Ve), grad(u_Ve))+(mu+ll)*inner(div(u_Ve), div(u_Ve)))*dx-inner(Constant((0,0)),u_Ve)*dx+inner(Constant((0.1,0)),u_Ve)*dd(1)
    err = assemble(error)
    print "RESIDUAL = ", err

    if u_e is not None:
        l2 = errornorm(u_e, uu, norm_type='L2', degree_rise=3)
        ene = errornorm(u_e, uu, norm_type='H10', degree_rise=3)
        print "L2 error = ", l2, " Energieerror = ", ene
    return uu




dp = 0.05
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

p2 = np.array([[0,0],[1,0],[0,1],[1,0.5], [0.5,1], [0.5,0.5]])

fd = lambda p: dm.dpoly(p, pv)

fd2 = lambda p: np.maximum(dm.drectangle0(p, 0,1, 0,1), -dm.drectangle0(p, 0.5,1.1, 0.5,1.1))

meshes = []

meshes.append([TriangularMesh.create_uniform_mesh(fd, [0, 1, 0, 1],0.01), "ref_mesh"])
meshes.append([TriangularMesh.create_bad_uniform_mesh(fd, [0,1,0,1], 0.05), "bad uniform"])
meshes.append([TriangularMesh.create_uniform_mesh(fd, [0, 1, 0, 1],0.04), "uniform"])
meshes.append([TriangularMesh.create_elliptic_mesh(20, Expression("4*x[0]+x[1]"), Expression("x[0]+4*x[1]"), fd), "elliptic"])
meshes.append([TriangularMesh.create_quasi_random_mesh(fd, [0,1,0,1], 350, 2, 3, pv), "quasi-random"])
meshes.append([TriangularMesh.create_spring_based_mesh(fd, [0,1,0,1], 0.05, pv, dm.huniform), "uniform_distmesh"])

fh = lambda p: 7.5*(p[:,0]-0.5)**2+7.5*(p[:,1]-0.5)**2+1
meshes.append([TriangularMesh.create_spring_based_mesh(fd, [0,1,0,1], 0.03, pv, fh), "corner_distmesh"])

u_e = None

for mesh, name in meshes:
    print name, len(mesh.get_nodes())
    print "minimum angle " , compute_minimum_angle(mesh)
    print "aspect ratio " , compute_aspect_ratio(mesh)
    print "edge ratio " , compute_edge_ratio(mesh)
    print "skewness " , compute_skewness(mesh)
    if "distmesh" in name:
        mesh.save_to_file(name, eps=1e-6)
    else:
        mesh.save_to_file(name)

    if name == "ref_mesh":
        u_e = solve_problem(name, u_e, p=5)
    else:
        solve_problem(name, u_e)
    print "------------------------------------------------------------"
    print ""

