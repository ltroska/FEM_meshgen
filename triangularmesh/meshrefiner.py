import numpy as np
from utils import compute_diameter, length2
from face import Face

__author__ = 'lukas'

class MeshRefiner(object):

    def __delete_or_add_non_conforming_edge(self, edge, edgeList):
        if edge:
            if edge[0] >= 0:
                edgeList.append(edge)
            else:
                try:
                    edgeList.remove([-edge[i] for i in range(len(edge))])
                except: pass

    def __split_triangle(self, x, y, z, p, t, old_triangle_index):
        """split side xz in the middle
            edge contains triangle index of x,y,v
        """
        px = p[x]
        py = p[y]
        pz = p[z]

        xz = pz-px
        pv = px+0.5*xz

        l = np.where(np.all(p==pv, axis=1))[0]

        if not l:
            v = len(p)
            p = np.append(p, [pv], axis=0)
            non_conforming_edge = [x, z, v, len(t)-1]
        else:
            v = l[0]
            non_conforming_edge = [-x, -z, -v, -len(t)+1]

        parent_diameter = t[old_triangle_index].diameter
        parent_index = t[old_triangle_index].index

        t = np.delete(t, old_triangle_index, 0)
        t = np.append(t, [Face([x,y,v], compute_diameter(px, py, pv), parent_diameter, parent_index)], axis=0)
        t = np.append(t, [Face([y,z,v], compute_diameter(py, pz, pv), parent_diameter, parent_index)], axis=0)

        return p, t, non_conforming_edge

    def __divide_by_longest_edge(self, triangle, p, t):
        x, y, z = t[triangle]
        px = p[x]
        py = p[y]
        pz = p[z]

        xy = py-px
        xz = pz-px
        yz = pz-py

        lxy = length2(xy)
        lxz = length2(xz)
        lyz = length2(yz)

        if lxy >= lxz and lxy >= lyz:
            p, t, non_conforming_edge = self.__split_triangle(x, z, y, p, t, triangle)
        elif lxz >= lxy and lxz >= lyz:
            p, t, non_conforming_edge= self.__split_triangle(x, y, z, p, t, triangle)
        elif lyz >= lxy and lyz >= lxz:
            p, t, non_conforming_edge = self.__split_triangle(y, x, z, p, t, triangle)

        return p, t, non_conforming_edge

    def __fix_non_conforming_edge(self, edge, p, t):
        x, y, v, old_triangle_index = edge

        #look for triangles containing this edge
        triangleIndex = -1
        for i in range(len(t)):
            if x in t[i] and y in t[i]:
                triangleIndex = i
                break

        #no triangle found
        if triangleIndex == -1:
            return p, t, None

        #get triangle
        triangle = t[triangleIndex]

        for vertex in triangle:
            if vertex != x and vertex != y:
                z = vertex

        px, py, pz = p[x], p[y], p[z]

        xy = py-px
        xz = pz-px
        yz = pz-py

        lxy = length2(xy)
        lxz = length2(xz)
        lyz = length2(yz)

        if lxy >= lxz and lxy >= lyz:
            p, t, _ = self.__split_triangle(x, z, y, p, t, triangleIndex)
            non_conforming_edge = None
        elif lxz >= lxy and lxz >= lyz:
            p, t, non_conforming_edge = self.__split_triangle(x, y, z, p, t, triangleIndex)
            p, t, _ = self.__split_triangle(x, abs(non_conforming_edge[2]), y, p, t, abs(non_conforming_edge[3]))
        elif lyz >= lxy and lyz >= lxz:
            p, t, non_conforming_edge = self.__split_triangle(y, x, z, p, t, triangleIndex)
            p, t, _ = self.__split_triangle(x, abs(non_conforming_edge[2]), y, p, t, abs(non_conforming_edge[3]))

        return p, t, non_conforming_edge


    def refine_by_longest_edge(self, triangles, p, t):
        non_conforming_edges = []
        triangleNodes = []

        for triangle in triangles:
            triangleNodes.append(t[triangle])

        for tr in triangleNodes:
            for triangle in range(len(t)):
                if t[triangle] == tr:
                    break
            p, t, non_conforming_edge = self.__divide_by_longest_edge(triangle, p, t)

            self.__delete_or_add_non_conforming_edge(non_conforming_edge, non_conforming_edges)

        while non_conforming_edges:
            non_conforming_edge = non_conforming_edges.pop(0)
            p, t, non_conforming_edge2 = self.__fix_non_conforming_edge(non_conforming_edge, p, t)

            self.__delete_or_add_non_conforming_edge(non_conforming_edge2, non_conforming_edges)

        return p, t
