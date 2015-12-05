__author__ = 'lukas'

from triangularmesh.mesh import *
from dolfin import *
import math


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
fd2= lambda p: np.maximum(dm.drectangle0(p, 0,1, 0,1), -dm.drectangle0(p, 0.5,1, 0.5,1))

P = Expression("4*x[0]+x[1]")
Q = Expression("x[0]+4*x[1]")

mesh = TriangularMesh.create_elliptic_mesh(15, P, Q, fd2)
fmesh = mesh.to_dolfin_mesh()

file = File("meshes/elliptic_example.pvd")
file << fmesh

plot(fmesh)
interactive()