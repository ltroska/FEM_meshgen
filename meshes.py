from triangularmesh.mesh import *
import math



meshes = []

h = .25
h2 = .2

for i in range(5):
    dp = h
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
    fh = lambda p: 7.5*(p[:,0]-0.5)**2+7.5*(p[:,1]-0.5)**2+1

	
	TriangularMesh.create_quasi_random_mesh(fd, [0,1,0,1], int(ceil(1/(h*h))), 2, 3, pv)

    tmp = TriangularMesh.create_bad_uniform_mesh(fd, [0,1,0,1], h2)
    meshes.append([tmp, "bad uniform_"+str(len(tmp.get_nodes()))])

    tmp = TriangularMesh.create_uniform_mesh(fd, [0, 1, 0, 1],h)
    meshes.append([tmp, "uniform_"+str(len(tmp.get_nodes()))])

    tmp = TriangularMesh.create_elliptic_mesh(int(ceil(1/h)), Expression("4*x[0]+x[1]"), Expression("x[0]+4*x[1]"), fd)
    meshes.append([tmp, "elliptic_"+str(len(tmp.get_nodes()))])

    tmp = TriangularMesh.create_quasi_random_mesh(fd, [0,1,0,1], int(ceil(1/(h*h))), 2, 3, pv)
    meshes.append([tmp, "quasi_random_"+str(len(tmp.get_nodes()))])

    tmp = TriangularMesh.create_spring_based_mesh(fd, [0,1,0,1], h, p2, dm.huniform)
    meshes.append([tmp, "uniform_distmesh_"+str(len(tmp.get_nodes()))])

    tmp = TriangularMesh.create_spring_based_mesh(fd, [0,1,0,1], h, p2, fh)
    meshes.append([tmp, "corner_distmesh_"+str(len(tmp.get_nodes()))])

    h/=2
    h2/=2

for mesh, name in meshes:
    fmesh = mesh.to_dolfin_mesh()
    plot(fmesh)
    interactive()
    mesh.save_to_file("meshes/"+name)
