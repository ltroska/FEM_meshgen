from mesh import *
import math
from utils import length

def compute_edge_ratio(mesh):
    """
    :param mesh: a mesh
    :return: minimum, maximum and average edge ratio of the mesh
    """

    p = mesh.get_nodes()
    t = mesh.get_faces()

    maximum = 0
    avg = 0

    #big number !
    minimum = 1e10

    for triangle in t:
        a, b, c = triangle
        pa, pb, pc = p[a], p[b], p[c]

        ab = length(pa-pb)
        ac = length(pa-pc)
        bc = length(pb-pc)


        local_min = min(min(ab,ac), bc)
        local_max = max(max(ab,ac), bc)


        edge_ratio = local_max/local_min

        minimum = minimum if minimum <= edge_ratio else edge_ratio
        maximum = maximum if maximum >= edge_ratio else edge_ratio

        avg += edge_ratio

    return minimum, maximum, avg/len(t)


def compute_aspect_ratio(mesh):
    """
    :param mesh: a mesh
    :return: minimum, maximum and average aspect ratio of the mesh
    """

    p = mesh.get_nodes()
    t = mesh.get_faces()

    maximum = 0
    avg = 0

    #big number !
    minimum = 1e10

    for triangle in t:
        a, b, c = triangle
        pa, pb, pc = p[a], p[b], p[c]

        ab = length(pa-pb)
        ac = length(pa-pc)
        bc = length(pb-pc)


        half_circumference = (ab+ac+bc)/2.
        area = math.sqrt(half_circumference*(half_circumference-ab)*(half_circumference-ac)*(half_circumference-bc))

        if area == 0:
            print pa, pb, pc, triangle

        inradius = area/half_circumference

        local_max = max(max(ab,ac), bc)


        aspect_ratio = local_max/(2*math.sqrt(3)*inradius)

        minimum = minimum if minimum <= aspect_ratio else aspect_ratio
        maximum = maximum if maximum >= aspect_ratio else aspect_ratio

        avg += aspect_ratio

    return minimum, maximum, avg/len(t)

def compute_minimum_angle(mesh):
    """
    :param mesh: a mesh
    :return: minimum, maximum and average minimum angle of the mesh
    """

    p = mesh.get_nodes()
    t = mesh.get_faces()

    maximum = 0
    avg = 0

    minimum = 91

    for triangle in t:
        a, b, c = triangle
        pa, pb, pc = p[a], p[b], p[c]

        ab = pb-pa
        ab = ab/sqrt(ab.dot(ab))
        ac = pc-pa
        ac = ac/sqrt(ac.dot(ac))
        bc = pc-pb
        bc = bc/sqrt(bc.dot(bc))

        angle_a = math.acos(ab.dot(ac))
        angle_b = math.acos((-ab).dot(bc))
        angle_c = math.acos((-ac).dot(-bc))

        local_min = min(min(angle_a, angle_b), angle_c) *180/math.pi

        minimum = minimum if minimum <= local_min else local_min
        maximum = maximum if maximum >= local_min else local_min

        avg += local_min

    return minimum, maximum, avg/len(t)

def compute_skewness(mesh):
    """
    :param mesh: a mesh
    :return: minimum, maximum and average skewness of the mesh
    """

    p = mesh.get_nodes()
    t = mesh.get_faces()

    maximum = 0
    avg = 0

    minimum = 1e10


    for i, triangle in enumerate(t):
        a, b, c = triangle
        pa, pb, pc = p[a], p[b], p[c]

        ab = length(pa-pb)
        ac = length(pa-pc)
        bc = length(pb-pc)

        half_circumference = (ab+ac+bc)/2.
        area = math.sqrt(half_circumference*(half_circumference-ab)*(half_circumference-ac)*(half_circumference-bc))

        circumscribed_circle_r = ab*ac*bc/(4.*area)

        a = math.sqrt(3)*circumscribed_circle_r
        optimal_area = math.sqrt(3)/4. * a*a

        skewness = (optimal_area-area)/optimal_area


        minimum = minimum if minimum <= skewness else skewness
        maximum = maximum if maximum >= skewness else skewness

        avg += skewness

    return minimum, maximum, avg/len(t)
