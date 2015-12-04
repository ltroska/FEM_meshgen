__author__ = 'lukas'

from meshrefiner import *
from dolfin import *
from utils import unique2d, closest_node, halton_sequence, two_nearest
from face import Face
import distmesh as dm
import scipy.spatial as spspatial
from math import log, pow, ceil


class TriangularMesh(object):
    def __init__(self):
        self.node_list = np.array([])
        self.face_list = np.array([])
        self.__meshrefiner = MeshRefiner()
        self.full_refines = 0
        self.fd = None
        self.dolfin_mesh = None

    def set_mesh(self, nodes, faces=None):
        """nodes in coordinates, faces in indices"""
        self.node_list = np.array(nodes.copy())
        if faces is not None:
            self.face_list = np.array(
                [Face(face, compute_diameter(nodes[face[0]], nodes[face[1]], nodes[face[2]]), -1) for face in faces])

    def add_node(self, node):
        if not node in self.node_list:
            np.append(self.node_list, node, axis=0)

    def add_face(self, node1, node2, node3):
        if all([set([node1, node2, node3]) != set(face) for face in self.face_list]):
            np.append(self.face_list, Face([node1, node2, node3]), axis=0)

    def remove_face(self, face):
        try:
            self.face_list.remove(face)
        except:
            pass

    def refine_all_with_bigger_error(self, errors, limit, refine_method='longest edge'):
        """refines all elements that have bigger error than limit"""
        triangles_to_refine = []

        for i, triangle in enumerate(self.face_list):
            if errors[i] > limit:
                triangles_to_refine.append(i)

        if 'longest' in refine_method:
            self.refine_by_longest_edge(triangles_to_refine)
        else:
            self.refine_by_delaunay(triangles_to_refine)

    def refine_all(self, refine_method='longest edge'):
        """refines all elements"""
        if 'longest' in refine_method:
            self.refine_by_longest_edge(range(len(self.face_list)))
        else:
            self.refine_by_delaunay(range(len(self.face_list)))
        self.full_refines += 1

    def refine_by_longest_edge(self, faces):
        """refines given elements/faces by the longest edge method"""
        self.node_list, self.face_list = self.__meshrefiner.refine_by_longest_edge(faces, self.node_list,
                                                                                   self.face_list)

    def to_dolfin_mesh(self):
        """returns a mesh in dolfin/FEniCS usable format"""
        self.update_index()
        mesh = Mesh()
        editor = MeshEditor()
        editor.open(mesh, 'triangle', 2, 2)
        editor.init_vertices(len(self.node_list))
        editor.init_cells(len(self.face_list))

        for idx in range(len(self.node_list)):
            editor.add_vertex(idx, self.node_list[idx][0], self.node_list[idx][1])

        for idx in range(len(self.face_list)):
            editor.add_cell(idx, self.face_list[idx][0], self.face_list[idx][1], self.face_list[idx][2])

        editor.close()

        return mesh

    def save_to_file(self, filename):
        """saves mesh to "filename".xml, boundaries are saved to "filename"_facets.xml"""
        file1 = File(filename + '.xml')
        file2 = File(filename + '_facets.xml')

        mesh = self.to_dolfin_mesh()

        class LeftBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] < DOLFIN_EPS

        class BottomBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and (x[1] < DOLFIN_EPS)

        class RightBoundaryVertical(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] > 1 - DOLFIN_EPS and x[1] < 0.5

        class RightBoundaryHorizontal(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] > 0.5 and abs(x[1]-0.5) <= DOLFIN_EPS

        class TopBoundaryHorizontal(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[1] > 1- DOLFIN_EPS and x[0] < 0.5

        class TopBoundaryVertical(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[1] >= 0.5 - DOLFIN_EPS and x[0] >= 0.5 - DOLFIN_EPS

        # compute topology by calling mesh.topology()

        mesh.topology()
        sub_domains = MeshFunction("size_t", mesh, 1)

        sub_domains.set_all(0)
        leftboundary = LeftBoundary()
        leftboundary.mark(sub_domains, 1)
        bottomboundary = BottomBoundary()
        bottomboundary.mark(sub_domains, 2)
        rightboundaryvert = RightBoundaryVertical()
        rightboundaryvert.mark(sub_domains, 3)
        rightboundaryhoriz = RightBoundaryHorizontal()
        rightboundaryhoriz.mark(sub_domains, 4)
        topboundaryhoriz =TopBoundaryHorizontal()
        topboundaryhoriz.mark(sub_domains, 5)
        topboundaryvert =TopBoundaryVertical()
        topboundaryvert.mark(sub_domains, 6)

        file1 << mesh
        file2 << sub_domains

    def adaptive_refine(self, errors_prev, errors_now, ratio=0.2, acceptable_error=1e-5):
        """adaptively refines the mesh using the Babuska-Rheinboldt strategy"""
        triangles_to_refine = []
        num_triangles_needed = int(ceil(ratio * len(self.face_list)))

        for i, triangle in enumerate(self.face_list):
            if triangle.parent_index == -1:
                continue
            parent_error = errors_prev[triangle.parent_index]
            parent_diameter = triangle.parent_diameter

            child_error = errors_now[i]
            child_diameter = triangle.diameter

            if child_error < acceptable_error:
                continue

            if abs(parent_diameter - child_diameter) < 1e-6:
                child_diameter = parent_diameter / sqrt(2)
            try:
                lbda = (log(child_error) - log(parent_error)) / (log(child_diameter) - log(parent_diameter))
                c = parent_error / (pow(parent_diameter, lbda))
                if c * pow(child_diameter / 2, lbda) < child_error:
                    triangles_to_refine.append(c * pow(child_diameter / 2, lbda))
            except:
                pass

        triangles_to_refine = np.array(triangles_to_refine)
        triangles_to_refine = triangles_to_refine[triangles_to_refine < 2]
        largest_estimated_error = max(triangles_to_refine)

        triangles_to_refine = []

        for i, triangle in enumerate(self.face_list):
            if errors_now[i] > largest_estimated_error:
                triangles_to_refine.append(i)

        if len(triangles_to_refine) <= num_triangles_needed:
            self.dolfin_mesh = self.to_dolfin_mesh()

        try:
            while (len(triangles_to_refine) <= num_triangles_needed):
                for tr in triangles_to_refine:
                    triangle = self.face_list[tr]
                    for node in triangle:
                        for i, new_triangle in enumerate(self.face_list):
                            if node in new_triangle and i not in triangles_to_refine:
                                triangles_to_refine.append(i)
                                if len(triangles_to_refine) >= num_triangles_needed:
                                    # get out of loop
                                    raise Exception
        except:
            pass

        print "refining ", len(triangles_to_refine), " of ", len(self.face_list), " edges"

        self.refine_by_longest_edge(triangles_to_refine)

    def refine_by_delaunay(self, triangles):
        """refines given elements by the delaunay refinement method"""
        p = self.node_list
        t = self.face_list

        for tr in triangles:
            pmid = (p[t[tr][0]] + p[t[tr][1]] + p[t[tr][2]]) / 3
            p = np.append(p, [pmid], axis=0)
        t = spspatial.Delaunay(p).vertices
        # remove triangles where centroid outside boundary (good enough...)
        pmid = p[t].sum(1) / 3
        t = t[self.fd(pmid) < 1e-15]

        t = [Face(t[i]) for i in range(len(t))]
        self.node_list = p
        self.face_list = t

    def update_index(self):
        for i in range(len(self.face_list)):
            self.face_list[i].index = i

    def get_diameters(self):
        """returns the diameter of all the faces"""
        return [face.diameter for face in self.face_list]

    def get_faces(self):
        self.update_index()
        return np.array([face.nodes for face in self.face_list])

    def get_nodes(self):
        return self.node_list

    @staticmethod
    def create_uniform_mesh(distance_function, bounding_box, edge_length, fixed_points=None):
        """
        distance_function > 0 outside, <0 inside, =0 boundary
        bounding box = [xmin, xmax, ymin, ymax, ...]
        edge_length = desired edge length of triangle
        fixed points = clear
        """
        geps = .001 * edge_length

        xmin, xmax, ymin, ymax = bounding_box
        mesh = TriangularMesh()

        # make grid
        x, y = np.mgrid[xmin:(xmax + edge_length):edge_length,
               ymin:ymax + edge_length:edge_length]

        # shift even rows
        x[:, 1::2] += edge_length / 2

        p = np.vstack((x.flat, y.flat)).T

        # remove points outside boundary
        p = p[distance_function(p) < -geps]

        # add fixed points
        if fixed_points is not None:
            fixed_points = np.array(fixed_points, dtype='d')
            p = np.vstack((fixed_points, p))
            p = unique2d(p)

        # triangulation
        t = spspatial.Delaunay(p).vertices

        # remove triangles where centroid outside boundary (good enough...)
        pmid = p[t].sum(1) / 3
        t = t[distance_function(pmid) < geps]

        mesh.set_mesh(p, t)
        mesh.fd = distance_function
        return mesh

    @staticmethod
    def create_equilateral_uniform_mesh(distance_function, bounding_box, edge_length, fixed_points=None):
        """
        distance_function > 0 outside, <0 inside, =0 boundary
        bounding box = [xmin, xmax, ymin, ymax, ...]
        edge_length = desired edge length of triangle
        fixed points = clear
        """

        geps = .001 * edge_length

        xmin, xmax, ymin, ymax = bounding_box
        mesh = TriangularMesh()

        # make grid
        x, y = np.mgrid[xmin:(xmax + edge_length):edge_length,
               ymin:(ymax + edge_length * np.sqrt(3) / 2):edge_length * np.sqrt(3) / 2]

        # shift even rows
        x[:, 1::2] += edge_length / 2

        p = np.vstack((x.flat, y.flat)).T

        for i in range(int(ceil((ymax - ymin) / (edge_length * np.sqrt(3) / 2)))):
            p = np.append(p, [[xmin, i * np.sqrt(3) / 2 * edge_length]], axis=0)
            p = np.append(p, [[xmax, i * np.sqrt(3) / 2 * edge_length]], axis=0)
            p = np.append(p, [[0.5, i * np.sqrt(3) / 2 * edge_length]], axis=0)

        for i in range(int(ceil((xmax - xmin) / edge_length))):
            p = np.append(p, [[xmin + 0.5 + i * edge_length, ymin + 0.5]], axis=0)

        # remove points outside boundary
        p = p[distance_function(p) < -geps]

        # add fixed points
        if fixed_points is not None:
            fixed_points = np.array(fixed_points, dtype='d')
            p = np.vstack((fixed_points, p))
            p = unique2d(p)

        # triangulation
        t = spspatial.Delaunay(p).vertices

        # remove triangles where centroid outside boundary (good enough...)
        pmid = p[t].sum(1) / 3
        t = t[distance_function(pmid) < geps]
        # print distance_function(p[t].sum(1)/3)

        mesh.set_mesh(p, t)
        mesh.fd = distance_function
        return mesh


    @staticmethod
    def create_spring_based_mesh(distance_function, bounding_box, initial_edge_length, fixed_points=None,
                                 element_size_function=dm.huniform):
        """
        distance_function > 0 outside, <0 inside, =0 boundary
        bounding box = [xmin, xmax, ymin, ymax, ...]
        initial_edge_length = edge length to start the initial uniform mesh with
        fixed points = clear
        element size function = relative desired edge length e.g. h(x,y)=1+x => edges at x=1 are twices as long as edges
                                at x=0, default uniform edge length
        """
        mesh = TriangularMesh()
        bounding_box[1], bounding_box[2] = bounding_box[2], bounding_box[1]
        p, t = dm.distmesh2d(distance_function, element_size_function, initial_edge_length,
                             bounding_box, fixed_points)
        mesh.set_mesh(p, t)
        mesh.fd = distance_function
        return mesh

    @staticmethod
    def create_quasi_random_mesh(distance_function, bounding_box, number_of_points, base1, base2, fixed_points=None):
        """
        Creates a mesh using the Halton sequence and triangulating them by Delaunay.
        :param distance_function: signed distance function for the region
        :param bounding_box: bounding box for the region, list [xmin, xmax, ymin, ymax]
        :param number_of_points: number of desired interior points
        :param base1: base for sequence in x coordinate
        :param base2: base for sequence in y coordinate
        :param fixed_point: array of predefined fixed points
        :return: a mesh instance representing the quasi random mesh
        """

        geps = 1e-10
        xmin, xmax, ymin, ymax = bounding_box
        mesh = TriangularMesh()

        p = fixed_points if fixed_points is not None else np.array([])

        count = 0
        index = 0

        while count < number_of_points:
            pv = [halton_sequence(index, base1), halton_sequence(index, base2)]

            if distance_function([pv]) < - geps:
                p = np.append(p, [pv], axis=0)
                count += 1

            index += 1

        t = spspatial.Delaunay(p).vertices

        # remove triangles with centroid outside domain
        pmid = p[t].sum(1) / 3
        t = t[distance_function(pmid) < geps]

        mesh.set_mesh(p, t)
        mesh.fd = distance_function
        return mesh

    @staticmethod
    def create_random_mesh(distance_function, bounding_box, number_of_points, dist=0, fixed_points=None):
        """
        distance_function > 0 outside, <0 inside, =0 boundary
        bounding box = [xmin, xmax, ymin, ymax, ...]
        number of points = clear
        dist = minimum distance for each new node to keep from all other points, default 0
        fixed points = clear
        if finding a new point fails more than 50 times, dist is set to dist=dist/sqrt(2)
        """

        geps = 1e-10
        n = number_of_points
        xmin, xmax, ymin, ymax = bounding_box
        mesh = TriangularMesh()

        p = fixed_points
        i = 0
        failures = 0

        while i < n:
            pv = [np.random.uniform(xmin, xmax), np.random.uniform(ymin, ymax)]

            if (dist == 0 or closest_node(pv, p) >= dist) and distance_function([pv]) < -geps:
                p = np.append(p, [pv], axis=0)
                i += 1
            else:
                failures += 1

            if failures > 50:
                failures = 0
                dist = dist / sqrt(2)

        # remove duplicates
        p = unique2d(p)

        # triangulation
        t = spspatial.Delaunay(p).vertices

        # remove triangles where centroid outside boundary (good enough...)
        pmid = p[t].sum(1) / 3
        t = t[distance_function(pmid) < geps]

        mesh.set_mesh(p, t)
        mesh.fd = distance_function
        return mesh

    @staticmethod
    def create_elliptic_mesh(number_points, P, Q, distance_function):
        """
        Note: might generate crappy meshes, not really optimized algorithm.
        :param number_points: number of points in each dimension => total of number_points**2 points
        :param P: Expression() for step size control in x direction
        :param Q: Expression() for step size control in y direction
        :return: a mesh generated by the elliptic pde using P and Q
        """

        # Create mesh and define function space
        mesh = UnitSquareMesh(number_points-1, number_points-1)
        V = FunctionSpace(mesh, 'Lagrange',1)

        def left_boundary(x, on_boundary):
            return on_boundary and x[0] < DOLFIN_EPS

        def bottom_boundary(x, on_boundary):
            return on_boundary and (x[1] < DOLFIN_EPS)

        def right_boundary1(x, on_boundary):
            return on_boundary and x[0] > 1 - DOLFIN_EPS and x[1] <= 0.5

        def right_boundary2(x, on_boundary):
            return on_boundary and x[0] > 1-DOLFIN_EPS and x[1] >= 0.5

        def top_boundary1(x, on_boundary):
            return on_boundary and x[1] > 1- DOLFIN_EPS and x[0] <= 0.5

        def top_boundary2(x, on_boundary):
            return on_boundary and x[1] > 1- DOLFIN_EPS and x[0] >= 0.5

        def getchi():

            u0_left_boundary = Constant(0)
            u0_bottom_boundary = Expression("x[0]")
            u0_right_boundary1 = Constant(1)
            u0_right_boundary2 =Expression("1.5-x[1]")
            u0_top_boundary1 = Expression("x[0]")
            u0_top_boundary2 = Constant(0.5)

            B1 = DirichletBC(V, u0_left_boundary, left_boundary)
            B2 = DirichletBC(V, u0_bottom_boundary, bottom_boundary)
            B3 = DirichletBC(V, u0_right_boundary1, right_boundary1)
            B4 = DirichletBC(V,u0_right_boundary2, right_boundary2)
            B5 = DirichletBC(V, u0_top_boundary1, top_boundary1)
            B6 = DirichletBC(V,  u0_top_boundary2, top_boundary2)

            BCs = [B1, B2, B3,B4, B5,  B6]

            # Define variational problem
            chi = TrialFunction(V)
            v = TestFunction(V)
            f = P
            a = -1*inner(nabla_grad(chi), nabla_grad(v))*dx
            L = f*v*dx

            # Compute solution
            chi = Function(V)
            solve(a == L, chi, BCs)

            # Plot solution and mesh
           # plot(chi)
            n = V.dim()
            d = mesh.geometry().dim()

            dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
            dof_coordinates.resize((n, d))
            dof_x = dof_coordinates[:, 0]
            dof_y = dof_coordinates[:, 1]

            return np.array(zip(dof_x, dof_y)), chi.vector().array()
            # Hold plot

        def getnu():

            u0_left_boundary = Expression("x[1]")
            u0_bottom_boundary = Constant(0)
            u0_right_boundary1 = Expression("x[1]")
            u0_right_boundary2 = Constant(0.5)
            u0_top_boundary1 = Constant(1)
            u0_top_boundary2 = Expression("1.5-x[0]")

            B1 = DirichletBC(V, u0_left_boundary, left_boundary)
            B2 = DirichletBC(V, u0_bottom_boundary, bottom_boundary)
            B3 = DirichletBC(V, u0_right_boundary1, right_boundary1)
            B4 = DirichletBC(V,u0_right_boundary2, right_boundary2)
            B5 = DirichletBC(V, u0_top_boundary1, top_boundary1)
            B6 = DirichletBC(V,  u0_top_boundary2, top_boundary2)

            BCs = [B1, B2, B3,B4, B5,  B6]

            # Define variational problem
            nu = TrialFunction(V)
            v = TestFunction(V)
            f = Q
            a = -1*inner(nabla_grad(nu), nabla_grad(v))*dx
            L = f*v*dx

            # Compute solution
            nu = Function(V)
            solve(a == L, nu, BCs)

            n = V.dim()
            d = mesh.geometry().dim()

            dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
            dof_coordinates.resize((n, d))
            dof_x = dof_coordinates[:, 0]
            dof_y = dof_coordinates[:, 1]

            return np.array(zip(dof_x, dof_y)), nu.vector().array()
            # Hold plot

        chixy, x = getchi()
        nuxy, y = getnu()

        nx = number_points

        mesh = TriangularMesh()

        p = np.array(zip(x, y))

        def compare(item1, item2):
            coords1 = item1[0]
            coords2 = item2[0]

            if coords1[0] < coords2[0]:
                return -1
            elif coords1[0] > coords2[0]:
                return 1
            elif coords1[0] == coords2[0]:
                if coords1[1] < coords2[1]:
                    return -1
                elif coords1[1] > coords2[1]:
                    return 1
                elif coords1[1] == coords2[1]:
                    return 0


        geps = 1e-6

        xy , p = zip(*sorted(zip(chixy, p), cmp=compare))

        p = np.array(p)

        t = np.array([[0,0,0]])

        for i in range(nx-1):
            for j in range(nx-1):
                t = np.append(t, [[i*nx+j, i*nx+j+1, (i+1)*nx+j+1]], axis=0)
                t = np.append(t, [[i*nx+j, (i+1)*nx+j, (i+1)*nx+j+1]], axis=0)


        t = np.delete(t, [0,0,0], 0)


        hanging_nodes = []

        for i, pt in enumerate(p):
            if distance_function(np.array([pt])) > geps:
                count = 0
                for j, tr in enumerate(t):
                    if i in tr:
                        a, b, c = tr
                        if not (a == pt).all():
                            hanging_nodes.append(a)
                        if not (b == pt).all():
                            hanging_nodes.append(b)
                        if not (c == pt).all():
                            hanging_nodes.append(c)

                        t = np.delete(t, j-count, 0)
                        count += 1

        hanging_nodes = set(hanging_nodes)
        iter_nodes = hanging_nodes.copy()

        for node in iter_nodes:
            intr = False
            for tr in t:
                if node in tr:
                    intr = True

            if not intr:
                hanging_nodes.remove(node)

        pnew = [0.5, 0.5]
        p = np.append(p, [pnew], axis=0)
        pnew = len(p)-1


        iter_nodes = hanging_nodes.copy()

        for node in iter_nodes:
            if node in hanging_nodes:
                for tr in t:
                    if node in tr:
                        a, b, c = tr
                        if a != node and a in hanging_nodes:
                            t = np.append(t, [[pnew, node, a]], axis=0)
                            #hanging_nodes.remove(a)
                            hanging_nodes.remove(node)
                            break

                        if b != node and b in hanging_nodes:
                            t = np.append(t, [[pnew, node, b]], axis=0)
                            #hanging_nodes.remove(b)
                            hanging_nodes.remove(node)
                            break

                        if c != node and c in hanging_nodes:
                            t = np.append(t, [[pnew, node, c]], axis=0)
                            #hanging_nodes.remove(c)
                            hanging_nodes.remove(node)
                            break

        hanging_nodes = []
        for i in range(len(p)):
            found = False

            for tr in t:
                if i in t:
                    found = True
                    break

            if not found:
                hanging_nodes.append(i)

        for idx in range(len(hanging_nodes)):
            node = hanging_nodes[idx]
            p = np.delete(p, node, 0)
            for tr in t:
                if tr[0] > node:
                    tr[0]-=1
                if tr[1] > node:
                    tr[1]-=1
                if tr[2] > node:
                    tr[2] -=1

            for i in range(len(hanging_nodes)):
                if hanging_nodes[i] > node:
                    hanging_nodes[i] -= 1


        if False:
            p = np.append(p , [[0.5, 1]], axis=0)
            p = np.append(p, [[1,0.5]], axis=0)

            a, b = two_nearest([0.5, 1], p)
            t = np.append(t, [[a, b, len(p)-2]], axis=0)

            a, b= two_nearest([1,0.5], p)
            t = np.append(t, [[a, b, len(p)-1]], axis=0)

        # remove triangles with centroid outside domain
        pmid = p[t].sum(1) / 3
        t = t[distance_function(pmid) < geps]
        mesh.fd = distance_function

        mesh.set_mesh(p, t)

        return mesh
