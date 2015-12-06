import numpy as np
from dolfin import *

def length2(p):
    return np.sum(p**2)

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sqrt(np.sum((nodes - node)**2, axis=1))
    return np.min(dist_2)

def two_nearest(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sqrt(np.sum((nodes - node)**2, axis=1))
    this = np.argmin(dist_2)
    dist_2[this] += 100000
    first = np.argmin(dist_2)
    dist_2[first] += 100000
    second = np.argmin(dist_2)

    return first, second

def halton_sequence(index, base):
    result = 0
    f = 1.
    i = index
    while i > 0:
        f = f/base
        result = result + f * (i % base)
        i = i/base

    return result

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
    return compute_energy_errors_on_cells(u_e, u, Ve, mesh)


def compute_energy_errors_on_cells(u_e, u, Ve, mesh):
    energy_errors = []

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

    return energy_errors

def compute_l2_errors_on_cells(u_e, u, Ve, mesh):
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

        error = e_Ve**2*dx(1)
        l2_errors.append(sqrt(assemble(error)))

    return l2_errors

def compute_residual_errors_on_cells(u, f, Ve, mesh):
    residuals = []

    cell_domains = CellFunction('size_t', mesh)
    u_Ve = interpolate(u, Ve)

    for cell in cells(mesh):
        cell_domains.set_all(0)
        cell_domains[cell] = 1

        dx = Measure('dx')[cell_domains]

        error = (inner(nabla_grad(u_Ve), nabla_grad(u_Ve))-f)*dx(1, domain=mesh)
        err = assemble(error)
        residuals.append(abs(err))


    return residuals

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


def length(p):
  return np.sqrt(np.sum(p**2))



"""
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
                yield t"""
