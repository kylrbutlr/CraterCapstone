"""
This is the Spice module. Put 'ya hot spicy code here.
More documentation to come.
"""

# Import packages/modules here
import spiceypy as spice
import datetime as datetime
import enum
import numpy as np
import Util

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
    :type date: datetime
    :return: time in Ephemeris Time Format(float)
    """

    return spice.utc2et(date)


#######
# FOV #
#######


def is_ray_in_fov(instrument_id, ray_direction, rframe, abcorr, observer, et):
    """
    Determines if a specified ray is in FOV of a specified instrument at a given time
    :param instrument_id: STR Name or ID code of the instrument
    :param ray_direction: 3-Element Array of floats - Ray's direction vector
    :param rframe: STR body-fixed, body-centered frame for target body
    :param abcorr: STR Aberration correction flag
    :param observer: STR Name or ID code of the observer
    :param et: FLOAT time of the observation (seconds past J2000)
    :return: BOOLEAN visibility flag
    """
    return spice.fovray(instrument_id, ray_direction, rframe, abcorr, observer, et)


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


def surface_intercept(target: str, time: float, fixed_reference: str, correction: str,
                               observer_name: str, direction_reference: str, direction_vector: list,
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

    return spice.getfov(instrument_id, max_return)


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


def print_ver():
    """
    Check to see if Spiceypy is installed
    :return: String = "CSPICE_N0066"
    """
    print(spice.tkvrsn('TOOLKIT'))


if __name__ == '__main__':
    # Load the kernels
    load_kernel('../kernels/mk/ROS_OPS.TM')

    # Load the image metafile information
    time, size = Util.get_lbl_information(r'../67P/ros_cam1_20150311t043620.lbl')

    # Find the position of the ROS_NAVCAM at time image was taken, based around 67P comet
    positions, lightTimes = sp_kernel_position('ROS_NAVCAM-A', convert_utc_to_et(time), 'J2000', 'NONE', '67P/C-G')
    positions = np.asarray(positions).T

    # Orientation Matrix for each of them (may be incorrect)
    # Orientation of ROS_NAVCAM frame to 67P Frame
    orientation_matrix = position_transformation_matrix("ROS_NAVCAM-A", "67P/C-G_CK", convert_utc_to_et(time))

    # Multiply orientation_matrix by camera vector
    # something like this because I don't have the camera vectors

    # for i in range(len(camera)):
    #    camera[i] = np.matmul(camera[i], orientation_matrix)

    # Camera/Instrument information
    shape, frame, bsight, n, bounds = get_instrument_fov(spice.bodn2c("ROS_NAVCAM-A"))

    unload_all_kernels()

