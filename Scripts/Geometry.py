"""
This is the Geometry module, the powerhouse of the crater project. More documentation to come.
"""

# Import packages/modules here
from Spice import *
import spiceypy as spice
import numpy as np
import math
from spiceypy.utils.support_types import SpiceyError
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

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


def init_kernels():
    """
    Load all the kernels
    :return: void
    """
    load_kernel('../kernels/mk/ROS_OPS.TM')
    load_kernel('../kernels/cassini/naif0009.tls')
    load_kernel('../kernels/cassini/cas00084.tsc')
    load_kernel('../kernels/cassini/cpck05Mar2004.tpc')
    load_kernel('../kernels/cassini/020514_SE_SAT105.bsp')
    load_kernel('../kernels/cassini/981005_PLTEPH-DE405S.bsp')
    load_kernel('../kernels/cassini/030201AP_SK_SM546_T45.bsp')
    load_kernel('../kernels/cassini/04135_04171pc_psiv2.bc')
    load_kernel('../kernels/cassini/cas_v37.tf')
    load_kernel('../kernels/cassini/cas_iss_v09.ti')
    load_kernel('../kernels/cassini/phoebe_64q.bds')
    load_kernel('../kernels/dsk/ROS_CG_M001_OSPCLPS_N_V1.BDS')


def get_et_values(utc):
    """
    Get et values one and two, so we can vectorize str2et
    :param utc: array of two utc dates
    :return: etOne and etTwo
    """
    etOne = convert_utc_to_et(utc[0])
    etTwo = convert_utc_to_et(utc[1])
    return [etOne, etTwo]


def get_id_code(name):
    """
    Initial check of id code of kernel, so program can stop if not found
    :param name: name of kernel, ex 'Rosetta'
    :return: id code if found or SpiceyError if not found
    """
    try:
        nacid = spice.bodn2c(name)
    except SpiceyError:
        # Stop the program if the code was not found.
        print('Unable to locate the ID code for ' + name)
        raise
    return nacid


def get_view_directions_vector(nacid):
    """
    Returns view directions of body
    :param nacid: body id
    :return: view directions vector
    """
    frame, vector, number, bounds, bounds2 = find_fov(nacid)
    return vector


def find_length_for_vertex_array(nacid):
    """
    Returns length for a vertex array
    :param nacid: body id
    :return: length
    """
    frame, vector, number, bounds, bounds2 = find_fov(nacid)
    bounds_list = bounds2.tolist()
    bounds_list.append(number.tolist())
    length = np.linspace(bounds2[1][0], bounds2[1][1], 150)
    return length


def conv_lat_to_rec_in_vertex_array(length):
    """
    converts latitude coordinates to rectangular in an array of size length
    :param length: size of array
    :return: vertex array
    """
    vertex_array = []
    for i in length:
        for j in length:
            vertex_array.append(spice.latrec(10000000000, i, j))
    vertex_array = np.array(vertex_array)
    return vertex_array


def get_direction_arrays(vertex_array):
    """
    Given a array of coordinates, find the directions of each
    :param vertex_array: array of rectangular coordinates
    :return: direction_arrays: array of directional arrays
    """
    n_rays = vertex_array.shape[0]
    direction_arrays = []
    for n in range(n_rays):
        direction_arrays.append([0, 0, 1])
    direction_arrays = np.array(direction_arrays)
    return direction_arrays


def get_position_vectors(vertex_array, etOne, body, vector):
    """
    Function gets x,y,z position vectors
    :param vertex_array: array of rectangular coordinates
    :param etOne: et value 1
    :param body: body name, ex. 'Rosetta'
    :param vector: views_directions vector
    :return: x,y,z data arrays in an array
    """
    xdata = []
    ydata = []
    zdata = []

    for vec in vertex_array:
        # Call sincpt to determine coordinates of the intersection of this vector with the surface
        try:
            point, trgepc, srfvec, area = find_ray_surface_intercept('67P/C-G', etOne, '67P/C-G_CK', 'NONE', body,
                                                                     vector, vec)
            # Now, we have discovered a point of intersection.
            #TODO do these need to printed?
            print(' Position vector of surface intercept in the 67P/C-G_CK frame (km):')
            print(' X = {:16.3f}'.format(point[0]))
            print(' Y = {:16.3f}'.format(point[1]))
            print(' Z = {:16.3f}'.format(point[2]))
            xdata.append(point[0])
            ydata.append(point[1])
            zdata.append(point[2])
        except:
            raise
    return [xdata, ydata, zdata]


def Brute_Force(body, nacid, utc):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # get et values
    etVals = get_et_values(utc)
    etOne = etVals[0]

    # Initial check for kernel id code
    get_id_code(body)

    # get view directions vector
    vector = get_view_directions_vector(nacid)

    # find length for vector array
    length = find_length_for_vertex_array(nacid)

    # fill in vertex array
    vertex_array = conv_lat_to_rec_in_vertex_array(length)

    # get position vectors for ray surface intercept
    position_vectors = get_position_vectors(vertex_array, etOne, body, vector)

    # plot position vectors
    ax.scatter3D(position_vectors[0], position_vectors[1], position_vectors[2], c=position_vectors[2])
    plt.show()

    if frame == 'RECTANGLE':  # will make a more robust version later on
        """dist1 = distance_formula(bounds2[0], bounds2[1])
        dist2 = distance_formula(bounds2[2], bounds2[3])
        area = dist1*dist2
        print(area)"""


        # print(vertex_array)
        # print(direction_arrays)
        # print(n_rays)



def dskxsi():
    """
    @Author: Kyler
    Will fix this hot garbage later
    :return:
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ROOM = 4
    nacid = None

    get_et_values(['2004-06-20', '2005-12-01'])

    nacid = get_id_code('CASSINI_ISS_NAC')

    frame, vector, number, bounds, bounds2 = find_fov(nacid, ROOM)
    #
    # Set values of "method" string that specify use of
    # ellipsoidal and DSK (topographic) shape models.
    #
    # In this case, we can use the same methods for calls to both
    # spiceypy.sincpt and spiceypy.ilumin. Note that some SPICE
    # routines require different "method" inputs from those
    # shown here. See the API documentation of each routine
    # for details.
    #
    method = ['Ellipsoid', 'DSK/Unprioritized']
    length = np.linspace(bounds2[1][0], bounds2[1][1], 10)
    vertex_array = []
    r = 10000000000
    for i in length:
        for j in length:
            vertex_array.append([i, j, bounds2[0][2]])
    # spice.vminus()
    # spice.dskgd
    # tup = spice.kdata(0, )
    vertex_array = np.array(vertex_array)
    n_rays = vertex_array.shape[0]
    direction_arrays = []
    print(vertex_array)
    for n in range(n_rays):
        direction_arrays.append([1, 0, 0])
    direction_arrays = np.array(direction_arrays)
    vertex_array = vertex_array.tolist()
    direction_arrays = direction_arrays.tolist()

    bounds_list = bounds2.tolist()
    bounds_list.append(number.tolist())
    #
    # Call sincpt to determine coordinates of the
    # intersection of this vector with the surface
    # of Phoebe.
    #
    tup = None
    try:
        tup = spice.dskxv(False, 'PHOEBE', [], 140254384.185, 'IAU_PHOEBE', vertex_array, direction_arrays)
        print(tup)
        xdata = []
        ydata = []
        zdata = []

        for vec in tup[0]:
            xdata.append(vec[0])
            ydata.append(vec[1])
            zdata.append(vec[2])
        ax.scatter3D(xdata, ydata, zdata)
        plt.show()
        # Now, we have discovered a point of intersection.
        # Start by displaying the position vector in the
        # IAU_PHOEBE frame of the intersection.
    except:
        raise


def distance_formula(a, b):
    x1 = a[0]
    y1 = a[1]
    x2 = b[0]
    y2 = b[1]
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

# Testing This stuff
ROSETTA = 'Rosetta'
ROSETTA_UTC = ['2016-12-31', '2016-12-31']
ROSETTA_NACID = -226807

if __name__ == '__main__':

    init_kernels()

    Brute_Force(ROSETTA, ROSETTA_NACID, ROSETTA_UTC)
    #dskxsi()

