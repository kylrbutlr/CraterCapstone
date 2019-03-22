"""
This is the Spice module. Put 'ya hot spicy code here.
More documentation to come.
"""

# Import packages/modules here
import spiceypy as spice
import enum
import numpy as np
import Util
import Geometry
import boundary
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
from spiceypy.utils.support_types import SpiceyError


def get_id_code(name):
    """
    check and return id code of kernel
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


def convert_lat_to_rec_coord(radius, longitude, latitude):
    """
    Convert from latitudinal coordinates to rectangular coordinates.
    :param radius: Distance of a point from the origin
    :param longitude: Longitude of point in radians.
    :param latitude: Latitude of point in radians.
    :return: Rectangular coordinates of the point. 3-element array of floats
    """
    return spice.latrec(radius, longitude, latitude)


###########
# Kernels #
###########


class KernelType(enum.Enum):
    """
    Enum types for kernels
    ** means kernels we will most likely be working with

    for more info about each type refer to:
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/12_intro_to_kernels.pdf
    """
    SPK = 'spk'  # **Spacecraft and Planet Ephemeris
    PCK = 'pck'  # Planetary Constants
    IK = 'ik'  # **Instrument
    CK = 'ck'  # **Camera Matrix (Orientation)
    EK = 'ek'  # Events
    FK = 'fk'  # **Reference Frame Specification
    SCLK = 'sclk'  # **Spacecraft Clock Correlation
    LSK = 'lsk'  # **Leap Seconds (Time)
    MK = 'mk'  # **Meta-Kernel (for loading many kernels)
    DSK = 'dsk'  # Digital Shape Kernel
    DBK = 'dbk'  # Database Mechanism


def load_kernel(file_name: str):
    """
    Load single or multiple kernels (through meta-kernel file, .mk)

    :param file_name: str file path of kernel
    """
    spice.furnsh(file_name)


def unload_kernel(file_name: str):
    """
    Unload a single or multiple kernels (through meta-kernel file, .mk)
    :param file_name: str file path of kernel
    """
    spice.unload(file_name)


def unload_all_kernels():
    """
    Unloads all kernels
    """
    spice.kclear()


def kernels_total(ktype: KernelType, value: int):
    """
    Check how many kernels are loaded for testing purposes
    :param ktype: Kernel enum type
    :param value: amount of kernels that should be loaded
    :return: boolean if there is a correct amount of kernels loaded
    """
    return spice.ktotal(ktype.value) == value


########
# Time #
########

def convert_utc_to_et(date: str):
    """
    Converts standard UTC date into ET (Ephemeris Time)

    :param date: input time in UTC form
    :type date: YYYY-MM-DDTHH:MM:SS.SSS = Year-Month-DayTHour:Minute:Second
    :return: time in Ephemeris Time Format(float)
    """

    return spice.utc2et(date)


def get_et_one(utc):
    """
    Get et date of first utc date
    :param utc: array of two utc dates
    :return: et date of first utc date
    """
    return convert_utc_to_et(utc[0])


def get_et_two(utc):
    """
    Get et date of second utc date
    :param utc: array of two utc dates
    :return: et date of second utc date
    """
    return convert_utc_to_et(utc[1])

###############
# SPK Kernels #
###############


def sp_kernel_position(main_body: str, time: float, reference_frame: str, correction: str, observing_body: str):
    """
    Loaded from SPK file
    Returns the position of a body relative to another observing body

    Example: sp_kernel_position('Cassini', 595166468.1826606, 'J2000', 'NONE', 'SATURN BARYCENTER')

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html
    :param main_body: the name of target body
    :param time: Ephemeris Time of the observer's location
    :param reference_frame: string referring to the reference frame of the output position vector
    :param correction:  Aberration correction flag
    :param observing_body: Name of observing body
    :return: Position of target in a 3D vector,
             Light time between the objects as a double
    """

    return spice.spkpos(main_body, time, reference_frame, correction, observing_body)


def position_transformation_matrix(from_object: str, to_object: str, time: float):
    """
    Returns the transformation matrix of how the frame is moving

    Example: position_transformation_matrix ( 'IAU_EARTH', 'J2000', ET  )

    ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pxform.html

    :param from_object: name of FK of which we are coming from
    :param to_object: name of FK of which we are going to
    :param time: Ephemeris Time Epoch of the rotation matrix
    :return: ET is the amount the position has transformed (magnitude?),
             The 3x3 transformation matrix M, where (x,y,z)^T *M = (x2,y2,z2)^T
    """

    return spice.pxform(from_object, to_object, time)


##############
# DSK Kernel #
##############

def find_ray_surface_intercept_at_epoch(target: str, time: float, fixed_reference: str, correction: str,
                                        observer_name: str, direction_reference:str, direction_vector: list,
                                        method: str = 'DSK/UNPRIORITIZED'):
    """
    Finds where a ray intercepts a surface in Cartesian coordinates

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/sincpt.html

    :param target: name of the target body
    :param time: Ephemeris Time of the frame
    :param fixed_reference: name of the FK that the target body is fixed to
    :param correction: Aberration correction
    :param observer_name: the name of the observer taking the frame
    :param direction_reference: name of FK that is the reference frame of which the ray's direction is
    :param direction_vector: Ray direction vector that is from the observer
    :param method: computational method name such as 'DSK/UNPRIORITIZED[/SURFACES = <surface list>]' or 'ELLIPSOID'
    :return: The closest point,
             The distance between the observer and the point,
             The observer's position at ET to the corrected position of surface intersect point,
             Boolean for if the ray intersects
    """

    return spice.sincpt(method, target, time, fixed_reference, correction, observer_name, direction_reference,
                        direction_vector)


def get_instrument_fov(instrument_id: int, max_return: int = 10):
    """
    Return information about the instrument such as the shape, reference frame, direction of view, and the corner
    of the instruments. The shape can be a 'POLYGON', 'RECTANGLE', 'CIRCLE', or 'ELLIPSE'.

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/25_ik.pdf
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html

    :param instrument_id: NAIF ID of the instrument (negative six digit number, example: -226801)
    :param max_return: maximum number of vectors that you want to return in 'BOUNDS' (recommended: 4+)
    :return:
             Shape of the fov,
             Reference frame in which the field of view boundary vectors are defined,
             Vector representing the view direction of the instrument,
             The number of 'BOUNDS' returned,
             'BOUNDS' that points to the corners of the instrument
    """
    shape, reference_frame, views_direction_vector, number_of_bounds, bounds = spice.getfov(instrument_id, max_return)
    return [shape, reference_frame, views_direction_vector, number_of_bounds, bounds]


def get_shape_fov(fov):
    return fov[0]


def get_reference_frame_fov(fov):
    return fov[1]


def get_views_direction_vector_fov(fov):
    return fov[2]


def get_number_of_bounds_fov(fov):
    return fov[3]


def get_bounds_fov(fov):
    return fov[4]


def get_shape_information(shapeFilePath: str):
    """
    Loads the shape file (dsk) and returns the vertices of the shape and
    the vertices associated to each plate

    :param shapeFilePath: relative path to the dsk files
    :return:
            Array of the coordinates of the vertices of the shape,
            Array of the indices of vertices that make up a plate
    """

    id = spice.dasopr(shapeFilePath)
    handle = spice.dlabfs(id)

    num_vertices, num_plates = spice.dskz02(id, handle)

    vertices_array = spice.dskv02(id, spice.dlabfs(id), 1, num_vertices)
    plate_array = spice.dskp02(id, handle, 1, num_plates)

    return vertices_array, plate_array


def distance_formula(p0, p1):
    return ((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)**(1/2)


def print_ver():
    """
    Check to see if Spiceypy is installed
    :return: String = "CSPICE_N0066"
    """
    print(spice.tkvrsn('TOOLKIT'))

def get_footprint(kernel_path:str, lbl_path:str, observing_body:str, target_body:str, target_body_frame:str):
    load_kernel(kernel_path)
    time, size = Util.get_lbl_information(lbl_path)
    size = int(size ** (1 / 2))
    footprint, points = Geometry.brute_force(observing_body, get_id_code(observing_body),
                                             convert_utc_to_et(time), size, target_body, target_body_frame)


    unload_all_kernels()
    return footprint

if __name__ == '__main__':

    region = get_footprint('../kernels/mk/ROS_OPS.TM', r'../67P/ros_cam1_20150408t061457.lbl', "ROS_NAVCAM-A",
                           "67P/C-G", "67P/C-G_CK")

    # Figure initialization
    ax = plt3d.Axes3D(plt.figure())
    ax.dist = 5
    ax.azim = -140
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])

    w = get_shape_information('../kernels/dsk/ROS_CG_M001_OSPCLPS_N_V1.bds')

    # Vertex data
    verts = w[0]
    # Face data
    faces = w[1]

    # Some triangles are weird so max is the threshold for triangle distance
    totalMax = 0

    # Draw the triangles
    for i in range(2500):
        max = 0
        num = np.random.randint(0, len(faces))

        # Get the points
        p1 = verts[faces[num, 0]]
        p2 = verts[faces[num, 1]]
        p3 = verts[faces[num, 2]]

        # Initialize the triangle
        triangle = [p1, p2, p3]
        face = plt3d.art3d.Poly3DCollection([triangle])
        face.set_color("grey")
        face.set_edgecolor('grey')
        face.set_alpha(.5)

        # Check to see if it passes the threshold
        a = distance_formula(p1, p2)
        b = distance_formula(p2, p3)
        c = distance_formula(p1, p3)

        # Collect the highest of the three numbers
        if a > b:
            max = a
        else:
            max = b
        if max < c:
            max = c

        # Was for seeing what a proper threshold is
        if max > totalMax:
            totalMax = max

        if max < 1.5:
            ax.add_collection3d(face)

    # Was for seeing what a proper threshold is

    for i in range(len(region)):
        xs = region[i, 0, 0], region[i, 1, 0]
        ys = region[i, 0, 1], region[i, 1, 1]
        zs = region[i, 0, 2], region[i, 1, 2]
        line = plt3d.art3d.Line3D(xs, ys, zs)
        ax.add_line(line)

    plt.show()


    unload_all_kernels()


