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


def Brute_Force():
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ROOM = 4
    nacid = None

    load_kernel('../kernels/kernels/mk/ROS_OPS.TM')
    # load_kernel('../kernels/kernels/naif0009.tls')
    # load_kernel('../kernels/kernels/cas00084.tsc')
    # load_kernel('../kernels/kernels/cpck05Mar2004.tpc')
    # load_kernel('../kernels/kernels/020514_SE_SAT105.bsp')
    # load_kernel('../kernels/kernels/981005_PLTEPH-DE405S.bsp')
    # load_kernel('../kernels/kernels/030201AP_SK_SM546_T45.bsp')
    # load_kernel('../kernels/kernels/04135_04171pc_psiv2.bc')
    # load_kernel('../kernels/kernels/cas_v37.tf')
    # load_kernel('../kernels/kernels/cas_iss_v09.ti')
    # load_kernel('../kernels/kernels/phoebe_64q.bds')
    load_kernel('../kernels/kernels/dsk/ROS_CG_M001_OSPCLPS_N_V1.BDS')
    utc = ['2014-08-06', '2016-12-31']

    # get et values one and two, we could vectorize str2et
    etOne = convert_utc_to_et(utc[0])
    etTwo = convert_utc_to_et(utc[1])
    # kinfo = spice.kinfo('../kernels/kernels/mk/ROS_OPS.TM')

    # print(spice.namfrm(frame_name))
    # tup = find_ray_surface_intercept('Rosetta', 595166468.1826606, '67P/C-G_FIXED', 'NONE', vector, 'J2000', )
    try:
        nacid = spice.bodn2c('Rosetta')
    except SpiceyError:
        # Stop the program if the code was not found.
        #
        print('Unable to locate the ID code for '
              'CASSINI_ISS_NAC')
        raise

    vecnam = ['Boundary Corner 1',
              'Boundary Corner 2',
              'Boundary Corner 3',
              'Boundary Corner 4',
              'Cassini NAC Boresight']
    frame, vector, number, bounds, bounds2 = find_fov(-226807)
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
    try:
        phoeid = spice.bodn2c('67P/C-G')
    except:
        #
        # The ID code for PHOEBE is built-in to the library.
        # However, it is good programming practice to get
        # in the habit of handling exceptions that may
        # be thrown when a quantity is not found.
        #
        print('Unable to locate the body ID code '
              'for Phoebe.')
        raise
        #
        # Now perform the same set of calculations for each
        # vector listed in the BOUNDS array. Use both
        # ellipsoidal and detailed (DSK) shape models.
        #

    # frame_name = spice.frmnam(1000012)
    bounds_list = bounds2.tolist()
    bounds_list.append(number.tolist())
    length = np.linspace(bounds2[1][0], bounds2[1][1], 150)
    vertex_array = []

    for i in length:
        for j in length:
            vertex_array.append(spice.latrec(1000000000.0, i, j))
    # spice.vminus()
    # spice.dskgd
    # tup = spice.kdata(0, )
    vertex_array = np.array(vertex_array)
    n_rays = vertex_array.shape[0]
    direction_arrays = []
    for n in range(n_rays):
        direction_arrays.append([0, 0, 1])
    direction_arrays = np.array(direction_arrays)
    xdata = []
    ydata = []
    zdata = []

    for vec in vertex_array:
            #
            # Call sincpt to determine coordinates of the
            # intersection of this vector with the surface
            # of Phoebe.
            #
        # print('Vector: {:s}\n'.format(vecnam[i]))
        try:
            point, trgepc, srfvec, area = find_ray_surface_intercept('67P/C-G', etTwo, '67P/C-G_CK', 'NONE',
                                                                         'Rosetta', vector, vec)
            #
            # Now, we have discovered a point of intersection.
            # Start by displaying the position vector in the
            # IAU_PHOEBE frame of the intersection.
            #
            print(' Position vector of surface intercept in the IAU_PHOEBE frame (km):')
            print(' X = {:16.3f}'.format(point[0]))
            print(' Y = {:16.3f}'.format(point[1]))
            print(' Z = {:16.3f}'.format(point[2]))
            xdata.append(point[0])
            ydata.append(point[1])
            zdata.append(point[2])
        except:
            raise

    ax.scatter3D(xdata, ydata, zdata, c=zdata)
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

    # sp.load_kernel('../kernels/kernels/mk/ROS_OPS.TM')
    load_kernel('../kernels/kernels/naif0009.tls')
    load_kernel('../kernels/kernels/cas00084.tsc')
    load_kernel('../kernels/kernels/cpck05Mar2004.tpc')
    load_kernel('../kernels/kernels/020514_SE_SAT105.bsp')
    load_kernel('../kernels/kernels/981005_PLTEPH-DE405S.bsp')
    load_kernel('../kernels/kernels/030201AP_SK_SM546_T45.bsp')
    load_kernel('../kernels/kernels/04135_04171pc_psiv2.bc')
    load_kernel('../kernels/kernels/cas_v37.tf')
    load_kernel('../kernels/kernels/cas_iss_v09.ti')
    load_kernel('../kernels/kernels/phoebe_64q.bds')
    utc = ['2004-06-20', '2005-12-01']
    etOne = convert_utc_to_et(utc[0])
    etTwo = convert_utc_to_et(utc[1])
    # kinfo = spice.kinfo('../kernels/kernels/mk/ROS_OPS.TM')

    # print(spice.namfrm(frame_name))
    # tup = find_ray_surface_intercept('Rosetta', 595166468.1826606, '67P/C-G_Ck', 'NONE', vector, 'J2000', )
    try:
        nacid = spice.bodn2c('CASSINI_ISS_NAC')
    except SpiceyError:
        # Stop the program if the code was not found.
        #
        print('Unable to locate the ID code for '
              'CASSINI_ISS_NAC')
        raise

    vecnam = ['Boundary Corner 1',
              'Boundary Corner 2',
              'Boundary Corner 3',
              'Boundary Corner 4',
              'Cassini NAC Boresight']
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
    try:
        phoeid = spice.bodn2c('PHOEBE')
    except:
        #
        # The ID code for PHOEBE is built-in to the library.
        # However, it is good programming practice to get
        # in the habit of handling exceptions that may
        # be thrown when a quantity is not found.
        #
        print('Unable to locate the body ID code '
              'for Phoebe.')
        raise
        #
        # Now perform the same set of calculations for each
        # vector listed in the BOUNDS array. Use both
        # ellipsoidal and detailed (DSK) shape models.
        #

    # frame_name = spice.frmnam(1000012)
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


if __name__ == '__main__':
    Brute_Force()
    # dskxsi()

