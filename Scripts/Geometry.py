"""
This is the Geometry module, the powerhouse of the crater project. More documentation to come.
"""

# Import packages/modules here
from Spice import *
import numpy as np
import boundary
from spiceypy.utils.support_types import SpiceyError


def find_bound_extrema(fovBounds):
    """
    Returns the min and max values of X and Y values of the fov bounds
    :param fovBounds: the fov bounds of the instrument
    :return: minX: the minimum value of X values in fov bounds
                maxX: the maximum value of X values in fov bounds
                minY: the minimum value of Y values in fov bounds
                maxY: the maximum value of Y values in fov bounds
    """
    # Initialize values
    minX = fovBounds[0][0]
    maxX = fovBounds[0][0]
    minY = fovBounds[0][1]
    maxY = fovBounds[0][1]

    for i in range(len(fovBounds)-1):
        x_value = fovBounds[i + 1][0]
        y_value = fovBounds[i + 1][1]

        if x_value > maxX:
            maxX = x_value
        if x_value < minX:
            minX = x_value

        if y_value > maxY:
            maxY = y_value
        if y_value < minY:
            minY = y_value

    return minX, maxX, minY, maxY


def get_direction_arrays(fovBounds, num_of_samples):
    """
    Returns the direction vectors within the fovBounds
    :param fovBounds: array of the edges of the field of view
    :param num_of_samples: the number of samples squared
    :return: direction_arrays: array of directional arrays
    """
    direction_arrays = []

    minX, maxX, minY, maxY = find_bound_extrema(fovBounds)
    diffX = maxX - minX
    diffY = maxY - minY

    # Evenly distribute points across the fov bounds
    for i in range(num_of_samples):
        for j in range(num_of_samples):
            direction_arrays.append([minX + i * (diffX / num_of_samples), minY + j * (diffY / num_of_samples),
                                     fovBounds[0][2]])

    return direction_arrays


def calculate_position_vectors_of_ray_surface_intercept(body, nacid, utc, num_of_samples, target_body,
                                                        target_body_reference, correction='NONE'):
    """
    Returns an array of many points of intersection at a specific epoch
    :param body: name of asteroid body
    :param nacid: id code of body
    :param utc: epoch
    :param num_of_samples: number of samples squared scanned within the field of view
    :param target_body: the naif name of the target body
    :param target_body_reference: the reference frame of the target body
    :param correction: aberration correction
    :return: x,y,z data arrays in an array
    """
    data = []
    fov = get_instrument_fov(nacid)
    reference_frame = get_reference_frame_fov(fov)
    if get_shape_fov(fov) == 'RECTANGLE':
        direction_array = get_direction_arrays(get_bounds_fov(fov), num_of_samples)

    for vec in direction_array:
        try:
            point, trgepc, srfvec = find_ray_surface_intercept_at_epoch(target_body, utc, target_body_reference,
                                                                        correction, body, reference_frame, vec)
            data.append(point)
        except SpiceyError:
            pass

    return data


def brute_force(body, nacid, utc, num_of_samples, target_body, target_body_reference):
    """
    Returns the footprint of the target body and the ray intercept points
    :param body: The observing or spacecraft body name
    :param nacid: The nacid value of the observing body
    :param utc: The epoch
    :param num_of_samples: The number of samples squared
    :param target_body: the naif name of the target body
    :param target_body_reference: the reference frame name of the target body
    :return: position_vectors: the coordinates of the edges of the footprint
                points: the collection of ray intercept points
    """

    position_vectors = calculate_position_vectors_of_ray_surface_intercept(body, nacid, utc, num_of_samples,
                                                                           target_body, target_body_reference)
    # plot_ray_surface_intercept(position_vectors)
    points = np.asarray(position_vectors)
    position_vectors = boundary.boundary(points)
    # plot_ray_surface_intercept(position_vectors)
    return position_vectors


#####################################################################################################################


def init_kernels():
    """
    Load all the kernels
    :return: void
    """
    load_kernel('../kernels/mk/ROS_OPS.TM')
    load_kernel('../kernels/dsk/ROS_CG_M001_OSPCLPS_N_V1.BDS')
    """
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
    l
    """



def find_length_for_vertex_array(fov, num_of_samples):
    """
    Returns length for a vertex array
    :param fov: field of view
    :param num_of_samples: number of samples for linspace
    :return: length
    """
    return np.linspace(get_bounds_fov(fov)[1][0], get_bounds_fov(fov)[1][1], num_of_samples)


def conv_lat_to_rec_in_vertex_array(length):
    """
    converts latitude coordinates to rectangular in an array of size length
    :param length: size of array
    :return: vertex array
    """
    vertex_array = []
    for i in length:
        for j in length:
            vertex_array.append(convert_lat_to_rec_coord(10000000000, i, j))
    vertex_array = np.array(vertex_array)
    return vertex_array


def create_vertex_array(fov, num_of_samples, rec=None):
    """
    Create an array for the coordinates of the fov
    :param fov: field of view
    :param num_of_samples: number of samples for linspace
    :param rec: determine whether or not to convert lat to rec or just append bounds to vertex array
    :return: array of rectangular coordinates
    """
    length = find_length_for_vertex_array(fov, num_of_samples)
    if rec == REC:
        vertex_array = conv_lat_to_rec_in_vertex_array(length)
    else:
        vertex_array = []
        for i in length:
            for j in length:
                vertex_array.append([i, j, get_bounds_fov(fov)[0][2]])
        vertex_array = np.array(vertex_array)
    return vertex_array


def calculate_intercept_point_array(nacid, ROOM, NUM_SAMPLES, TARGET, EPOCH, FIXEDREF):
    fov = find_fov(nacid, ROOM)
    vertex_array = create_vertex_array(fov, NUM_SAMPLES)
    direction_arrays = get_direction_arrays(vertex_array)
    vertex_array = vertex_array.tolist()
    direction_arrays = direction_arrays.tolist()

    intercept_array = find_ray_surface_intercept_by_DSK_segments(False, TARGET, [], EPOCH, FIXEDREF, vertex_array, direction_arrays)
    print(intercept_array)
    intercept_array = intercept_array[0]
    return intercept_array


def plot_ray_surface_intercept(intercept_vectors):
    """
    Plots the position vectors of the ray surface intercepts
    :param intercept_vectors: position vectors or intercept points
    :return: figure plot
    """

    xdata = []
    ydata = []
    zdata = []

    for intercept in intercept_vectors:
        xdata.append(intercept[0])
        ydata.append(intercept[1])
        zdata.append(intercept[2])

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(xdata, ydata, zdata, c=zdata)
    plt.show()


def dskxsi(body):
    """
    @Author: Kyler
    Will fix this hot garbage later
    :return:
    """
    nacid = get_id_code(body)
    intercept_point_array = calculate_intercept_point_array(nacid, CASSINI_ROOM, CASSINI_NUM_SAMPLES, CASSINI_TARGET, CASSINI_EPOCH, CASSINI_FIXEDREF)
    plot_ray_surface_intercept(intercept_point_array)


def distance_formula(a, b):
    x1 = a[0]
    y1 = a[1]
    x2 = b[0]
    y2 = b[1]
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

# Testing This stuff
ROSETTA = 'Rosetta'
ROSETTA_UTC = ['2016-12-31', '2016-12-31']
ROSETTA_NUM_SAMPLES = 150
ROSETTA_NACID = -226807
Target_Body = '67P/C-G'
Target_Body_ref_frame = '67P/C-G_CK'
CASSINI = 'CASSINI_ISS_NAC'
CASSINI_UTC = ['2004-06-20', '2005-12-01']
CASSINI_NUM_SAMPLES = 10
CASSINI_EPOCH = 140254384.185
REC = 'REC'
CASSINI_ROOM = 4
CASSINI_TARGET = 'PHOEBE'
CASSINI_FIXEDREF = 'IAU_PHOEBE'

if __name__ == '__main__':

    init_kernels()

    Brute_Force(ROSETTA, ROSETTA_NACID, ROSETTA_UTC, ROSETTA_NUM_SAMPLES, Target_Body, Target_Body_ref_frame)
    # dskxsi(CASSINI)
