"""
This is the Spice module. Put 'ya hot spicy code here.
More documentation to come.
"""

# Import packages/modules here
import spiceypy as spice
import datetime as datetime
import enum
from spiceypy.utils.support_types import SpiceyError

spice.dskxsi
spice.dskxv
def compute_view_frustrum(instrument):
    """
    stub
    :param instrument: Data structure that models the instrument.
    :return:
    """
    pass


def spacecraft_position_and_orientation(spacecraft, point):
    """
    Look under SPK Kernel for now
    :param spacecraft: Data structure that models the spacecraft.
    :param point: point in time.
    :return:
    """
    pass


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

"""
Create Spice Wrappers Here:
"""

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
    SPK = 'spk'     # **Spacecraft and Planet Ephemeris
    PCK = 'pck'     # Planetary Constants
    IK = 'ik'       # **Instrument
    CK = 'ck'       # **Camera Matrix (Orientation)
    EK = 'ek'       # Events
    FK = 'fk'       # **Reference Frame Specification
    SCLK = 'sclk'   # **Spacecraft Clock Correlation
    LSK = 'lsk'     # **Leap Seconds (Time)
    MK = 'mk'       # **Meta-Kernel (for loading many kernels)
    DSK = 'dsk'     # Digital Shape Kernel
    DBK = 'dbk'     # Database Mechanism


def load_kernel(file_name : str):
    """
    Load single or multiple kernels (through meta-kernel file, .mk)

    :param file_name: str file path of kernel
    """
    spice.furnsh(file_name)


def unload_kernel(file_name : str):
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


def check_number_kernels(ktype : KernelType, value: int):
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
    return convert_utc_to_et(utc[0])


###############
# SPK Kernels #
###############


def find_body_position(main_body: str, time: float, reference_frame: str, correction: str, observing_body: str):
    """
    Loaded from SPK file
    Returns the position of a body relative to another observing body

    Example: body_position('Cassini', 595166468.1826606, 'J2000', 'NONE', 'SATURN BARYCENTER')

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


def find_frame_transformation(from_object: str, to_object:str, time: float):
    """
    Returns the transformation matrix of how the frame is moving

    Example: find_frame_transformation ( 'IAU_EARTH', 'J2000', ET  )

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


def find_ray_surface_intercept_by_DSK_segments(pri, target, srflst, et, fixref, vtxarr, dirarr):
    """
    Compute ray-surface intercepts for a set of rays, using data provided by multiple loaded DSK segments.
    :param pri: (bool) – Data prioritization flag.
    :param target: (str) – Target body name.
    :param srflst: (list of int) – Surface ID list.
    :param et: (float) – Epoch, expressed as seconds past J2000 TDB.
    :param fixref: (str) – Name of target body-fixed reference frame.
    :param vtxarr: (Nx3-Element Array of floats) – Array of vertices of rays.
    :param dirarr: (Nx3-Element Array of floats) – Array of direction vectors of rays.
    :return: Intercept point array and Found flag array. Tuple
    """
    return spice.dskxv(pri, target, srflst, et, fixref, vtxarr, dirarr)


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


def find_fov(instrument_id: int, max_return: int = 10):
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


def print_ver():
    """
    Check to see if Spiceypy is installed
    :return: String = "CSPICE_N0066"
    """
    print(spice.tkvrsn('TOOLKIT'))

#Testing garbage out
if __name__ == '__main__':
    print_ver()
    load_kernel('../Tests/rosetta/kernels/mk/ROS_OPS.TM')
    print('break')
    print(spice.ktotal(KernelType.SPK.value))
    print(convert_utc_to_et(datetime.date.today()))
    frame, vector, number, bounds, bounds2 = find_fov(-226807)

    print(bounds2)
