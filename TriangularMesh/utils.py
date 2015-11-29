import numpy as np
from dolfin import *

def length(p):
    return p[0]*p[0]+p[1]*p[1]

class Face(object):
    def __init__(self, node_list, diameter = 0, parent_diameter = 0, parent_index = -1):
        if type(node_list) is not int:
            self.nodes = node_list
            self.parent_diameter = parent_diameter
            self.diameter = diameter
        else:
            self.nodes = [node_list, diameter, parent_diameter]
            self.parent_diameter = 0
            self.diameter = 0

        self.index = -1
        self.parent_index = parent_index

    def __eq__(self, other):
        return set(self) == set(other)

    def __str__(self):
        return '[' + str(self.nodes[0]) + ',' + str(self.nodes[1]) + ',' + str(self.nodes[2]) + ']'

    def __repr__(self):
        return '[' + str(self.nodes[0]) + ',' + str(self.nodes[1]) + ',' + str(self.nodes[2]) + ']'

    def __getitem__(self, item):
        return self.nodes[item]

    def __setitem__(self, key, value):
        self.nodes[key] = value

    class __metaclass__(type):
        def __iter__(self):
            for t in self.nodes:
                yield t

class Node(object):
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)
        self.faces = []

    def __eq__(self, other):
        return (self.coordinates == other.coordinates).all()

    def __str__(self):
        return str(self.coordinates)

    def __repr__(self):
        return str(self.coordinates)

    def __getitem__(self, item):
        return self.coordinates[item]

    def __setitem__(self, key, value):
        self.coordinates[key] = value

    class __metaclass__(type):
        def __iter__(self):
            for t in self.coordinates:
                yield t


def compute_diameter(px, py, pz):
    return max(length(px-py), length(px-pz), length(py-pz))

def energy_errornorm(u_e, u, Ve):
    u_Ve = interpolate(u, Ve)
    u_e_Ve = interpolate(u_e, Ve)
    e_Ve = Function(Ve)
    e_Ve.vector()[:] = u_e_Ve.vector().array()-u_Ve.vector().array()

    error = inner(nabla_grad(e_Ve), nabla_grad(e_Ve))*dx

    return sqrt(assemble(error))


def estimateErrorWithExactSolution(u_e, u, Ve, mesh):
    return computeErrorsOnCells(u_e, u, Ve, mesh)


def computeErrorsOnCells(u_e, u, Ve, mesh):
    energy_errors = []
    l2_errors = []

    u_Ve = interpolate(u, Ve)
    u_e_Ve = interpolate(u_e, Ve)
    e_Ve = Function(Ve)
    e_Ve.vector()[:] = u_e_Ve.vector().array()-u_Ve.vector().array()

    cell_domains = CellFunction('size_t', mesh)

    for cell in cells(mesh):
        cell_domains.set_all(0)
        cell_domains[cell] = 1

        dx = Measure('dx')[cell_domains]

        error = inner(nabla_grad(e_Ve), nabla_grad(e_Ve))*dx(1)
        energy_errors.append(sqrt(assemble(error)))

        error = e_Ve**2*dx(1)
        l2_errors.append(sqrt(assemble(error)))

    return energy_errors, l2_errors


def computeSolution(u0, f, mesh, degree, u_e = None, degree_rise = 0):
    V = FunctionSpace(mesh, 'Lagrange', degree)

    def u0_boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, u0, u0_boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(nabla_grad(u), nabla_grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    if u_e:
        return u, u_e
    elif degree_rise != 0:
        V = FunctionSpace(mesh, 'Lagrange', degree+degree_rise)
        u2 = TrialFunction(V)
        v = TestFunction(V)
        a = inner(nabla_grad(u2), nabla_grad(v))*dx
        L = f*v*dx
        u2 = Function(V)
        solve(a == L, u2, bc)

        return u, u2

    return u, None


def unique2d(a):
    x, y = a.T
    b = x + y*1.0j
    idx = np.unique(b,return_index=True)[1]
    return a[idx]

#PLOTTEN DER DIST FCT
"""x = np.arange(-1, 2, .01)
y = np.arange(-1,2, .01)
x, y = np.meshgrid(x,y)

r = np.array([[x[i,j], y[i,j]] for i in range(x.shape[0]) for j in range(x.shape[0])])

dist = fd(r)

dist[dist < 0] = -1
dist[dist > 0] = 1

dist = dist.reshape(x.shape)

plt.pcolormesh(x,y,dist, cmap='RdBu', vmin = -2, vmax=2)

plt.colorbar()
plt.show()"""

__author__ = 'lukas'

def length(p):
  return np.sqrt(np.sum(p**2))
