"""
This is the Geometry module, the powerhouse of the crater project. More documentation to come.
"""

# Import packages/modules here
import Spice as sp
import spiceypy as spice
import numpy as np
import math

def edges_around_plates(plates):
    """
    stub
    :param plates: Data structure that contains plates of the model.
    :return:
    """
    pass


def generate_texture_mapping(model):
    """
    stub
    :param model: The model without texture mapping.
    :return:
    """
    pass


def camera_vertices(camera):
    """
    stub
    :param camera: Data structure that contains information about the camera.
    :return:
    """
    pass


def plate_frustrum():
    pass


def view_frustrum_vertices():
    pass


def load_shape_data(file):
    """
    stub
    :param file: file location of teh shape data.
    :return:
    """
    pass


def Brute_Force():
    sp.load_kernel('../Tests/rosetta/kernels/mk/ROS_OPS.TM')
    frame, vector, number, bounds, bounds2 = sp.find_fov(-226807)

    tup = sp.find_ray_surface_intercept('Rosetta', 595166468.1826606, '1000012', 'NONE', vector, 'J2000', list(number))

    if frame == 'RECTANGLE':  # will make a more robust version later on
        """dist1 = distance_formula(bounds2[0], bounds2[1])
        dist2 = distance_formula(bounds2[2], bounds2[3])
        area = dist1*dist2
        print(area)"""

        length = np.linspace(bounds2[1][0], bounds2[1][1], 10)
        vertex_array = []

        for i in length:
            for j in length:
                vertex_array.append([i, j, bounds2[0][2]])
        # spice.vminus()
        # spice.dskgd
        # tup = spice.kdata(0, )
        vertex_array = np.array(vertex_array)
        n_rays = vertex_array.shape[0]
        direction_arrays = []
        for n in range(n_rays):
            direction_arrays.append([0, 0, 1])
        direction_arrays = np.array(direction_arrays)
        print(vertex_array)
        print(direction_arrays)
        print(n_rays)



def distance_formula(a, b):
    x1 = a[0]
    y1 = a[1]
    x2 = b[0]
    y2 = b[1]
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

# Testing This stuff


if __name__ == '__main__':
    Brute_Force()

