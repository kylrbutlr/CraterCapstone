"""
This is the Spice module. Put 'ya hot spicy code here.
More documentation to come.
"""

# Import packages/modules here
import spiceypy as spice
import datetime as datetime
import enum
import numpy as np
import matplotlib.pyplot as plt

def compute_view_frustum(instrument):
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


def convert_utc_to_et(date: datetime):
    """
    Converts standard UTC date into ET (Ephemeris Time)

    :param date: input time in UTC form
    :type date: datetime
    :return: time in Ephemeris Time Format(float)
    """

    w = datetime.date.strftime(date, '%b %d, %Y')
    return spice.utc2et(w)


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


def position_transformation_matrix(from_object: str, to_object: str, time: float):
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


def print_ver():
    """
    Check to see if Spiceypy is installed
    :return: String = "CSPICE_N0066"
    """
    print(spice.tkvrsn('TOOLKIT'))


# Testing garbage out
if __name__ == '__main__':
    print_ver()
    load_kernel('kernels/mk/ROS_OPS.TM')
    print('break')
    print(spice.ktotal(KernelType.SPK.value))
    print(convert_utc_to_et(datetime.date.today()))
    frame, vector, number, bounds, bounds2 = get_instrument_fov(-226807)

    print(bounds2)


    # Messing around with positions of Rosetta mission

    step = 4000
    # we are going to get positions between these two dates
    utc = ['2014-08-06', '2016-12-31']

    # get et values one and two, we could vectorize str2et
    etOne = spice.str2et(utc[0])
    etTwo = spice.str2et(utc[1])
    print("ET One: {}, ET Two: {}".format(etOne, etTwo))

    # get times
    times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

    # check first few times:
    print(times[0:3])

    positions, lightTimes = spice.find_body_position('ROSETTA', times, 'J2000', 'NONE', 'EARTH_BARYCENTER')

    # Positions is a 3xN vector of XYZ positions
    print("Positions: ")
    print(positions[0])

    # Light times is a N vector of time
    print("Light Times: ")
    print(lightTimes[0])

    positions = np.asarray(positions).T # positions is a list, make it an ndarray for easier indexing
    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111)
    ax.plot(positions[0], positions[1], positions[2])
    plt.title('Test')
    plt.show()

    unload_all_kernels()

    # Cassini Example from https://media.readthedocs.org/pdf/spiceypy/master/spiceypy.pdf
    load_kernel('kernels/naif0009.tls')
    load_kernel('kernels/cas00084.tsc')
    load_kernel('kernels/cpck05Mar2004.tpc')
    load_kernel('kernels/020514_SE_SAT105.bsp')
    load_kernel('kernels/981005_PLTEPH-DE405S.bsp')
    load_kernel('kernels/030201AP_SK_SM546_T45.bsp')
    load_kernel('kernels/04135_04171pc_psiv2.bc')
    load_kernel('kernels/cas_v37.tf')
    load_kernel('kernels/cas_iss_v09.ti')
    print(spice.ktotal(KernelType.SPK.value))

    step = 4000
    # we are going to get positions between these two dates
    utc = ['2004-06-20', '2005-12-01']

    # get et values one and two, we could vectorize str2et
    etOne = spice.str2et(utc[0])
    etTwo = spice.str2et(utc[1])
    print("ET One: {}, ET Two: {}".format(etOne, etTwo))

    # get times
    times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

    # check first few times:
    print(times[0:3])

    positions, lightTimes = spice.find_body_position('Cassini', times, 'J2000', 'NONE', 'SATURN BARYCENTER')

    # Positions is a 3xN vector of XYZ positions
    print("Positions: ")
    print(positions[0])

    # Light times is a N vector of time
    print("Light Times: ")
    print(lightTimes[0])

    positions = np.asarray(positions).T # positions is a list, make it an ndarray for easier indexing
    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111)
    ax.plot(positions[0], positions[1], positions[2])
    plt.title('SpiceyPy Cassini Position Example from Jun 20, 2004 to Dec 1, 2005')
    plt.show()
    unload_all_kernels()

