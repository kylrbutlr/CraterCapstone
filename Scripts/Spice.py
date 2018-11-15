"""
This is the Spice module. Put 'ya hot spicy code here.
More documentation to come.
"""

# Import packages/modules here
import spiceypy as spice
import datetime as datetime
import enum

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

"""
Create Spice Wrappers Here:
"""

###########
# Kernels #
###########

class KernelType(enum.Enum):
    """
    Enum types for kernels

    for more info about each type refer to:
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/12_intro_to_kernels.pdf
    """
    SPK = 'spk'     #**Spacecraft and Planet Ephemeris
    PCK = 'pck'     #Planetary Constants
    IK = 'ik'       #**Instrument
    CK = 'ck'       #**Camera Matrix (Orientation)
    EK = 'ek'       #Events
    FK = 'fk'       #**Reference Frame Specification
    SCLK = 'sclk'   #**Spacecraft Clock Correlation
    LSK = 'lsk'     #**Leap Seconds (Time)
    MK = 'mk'       #**Meta-Kernel (for loading many kernels)
    DSK = 'dsk'     #Digital Shape Kernel
    DBK = 'dbk'     #Database Mechanism

def load_kernel(fileName : str):
    """
    Load single or multiple kernels (through meta-kernel file, .mk)

    :param fileName: str file path of kernel
    """
    spice.furnsh(fileName)

def unload_kernel(fileName : str):
    """
    Unload a single or multiple kernels (through meta-kernel file, .mk)
    :param fileName: str file path of kernel
    """
    spice.unload(fileName)

def unload_all_kernels():
    """
    Unloads all kernels
    """
    spice.kclear()

def check_number_kernels(kType : KernelType, value: int):
    """
    Check how many kernels are loaded for testing purposes
    :param kType: Kernel enum type
    :param value: amount of kernels that should be loaded
    :return: boolean if there is a correct amount of kernels loaded
    """
    return spice.ktotal(kType.value) == value

########
# Time #
########

def convert_UTC_to_ET(date: datetime):
    """
    Converts standard UTC date into ET (Ephemeris Time)

    :param date: input time in UTC form
    :type date: datetime
    :return: time in Ephemeris Time Format(float)
    """

    w = datetime.date.strftime(date,'%b %d, %Y')
    return spice.utc2et(w)

###############
# SPK Kernels #
###############

def find_body_position(mainBody: str, etTime: float, referenceFrame: str, correction: str, observingBody: str):
    """
    Loaded from SPK file
    Returns the position of a body relative to another observing body

    Example: body_position('Cassini', 595166468.1826606, 'J2000', 'NONE', 'SATURN BARYCENTER')

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html
    :param mainBody: the name of target body
    :param Ephemeris Time: time in ET format (look at convert_UTC_to_ET)
    :param referenceFrame: string referring to the reference frame of the output position vector
    :param correction:  Abberation correction flag
    :param observingBody: Name of observing body
    :return: Position of target in a 3D vector
    :return: light time between the objects as a double
    """

    return spice.spkpos(mainBody, etTime, referenceFrame, correction, observingBody)

def find_frame_transformation(fromObject: str, toObject:str, etTime: float):
    """
    Returns the transformation matrix of how the frame is moving

    Example: find_frame_transformation ( 'IAU_EARTH', 'J2000', ET  )

    ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pxform.html

    :param fromObject: name of FK of which we are coming from
    :param toObject: name of FK of which we are going to
    :param Ephemeris Time: Epoch of the rotation matrix
    :return: ET is the amount the position has transformed (magnitude?)
    :return: The 3x3 transformation matrix M, where (x,y,z)^T *M = (x2,y2,z2)^T
    """

    return spice.pxform(fromObject, toObject, etTime)

##############
# DSK Kernel #
##############

def find_ray_surface_intercept(method: str, target: str, etTime: float, fixedReference:str, correction: str,
                  observerName: str, directionReference:str, directionVector: list):
    """
    Finds where a ray intercepts a surface

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/sincpt.html

    :param method: computational method name such as 'ELLIPSOID' or 'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
    :param target: name of the target body
    :param Ephemeris Time: et time of the frame
    :param fixedReference: name of the FK that the target body is fixed to
    :param correction: Aberration correction
    :param observerName: the name of the observer taking the frame
    :param directionReference: name of FK that is the reference frame of which the ray's direction is
    :param directionVector: Ray direction vector that is from the observer
    :return: The closest point in Cartesion coordinates
    :return: The distance between the observer and the intercept point
    :return: The observer's position at ET to the corrected position of surface intersect point
    :return: Boolean value which expresses if the ray intersects
    """

    return spice.sincpt(method, target, etTime, fixedReference, correction, observerName, directionReference,
                        directionVector)

def find_fov(instrumentID: int, maxReturn: int):
    """
    Return information about the instrument such as the shape, reference frame, direction of view, and the corner
    of the instruments. The shape can be a 'POLYGON', 'RECTANGLE', 'CIRCLE', or 'ELLIPSE'.

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/25_ik.pdf
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html

    :param instrumentID: NAIF ID of the instrument (negative six digit number, example: -226801)
    :param maxReturn: maximum number of vectors that you want to return in 'BOUNDS' (recommended: 4+)
    :return: shape of the fov
    :return: reference frame in which the field of view boundary vectors are defined
    :return: Vector representing the view direction of the instrument
    :return: the number of bounds returned
    :return: 'BOUNDS' that points to the corners of the instrument
    """

    return spice.getfov(instrumentID, maxReturn)


def print_ver():
    """
    Check to see if Spiceypy is installed
    :return: String = "CSPICE_N0066"
    """
    print(spice.tkvrsn('TOOLKIT'))

#Testing garbage out
if __name__ == '__main__':
    print_ver()
    load_kernel('kernels/mk/ROS_OPS.TM')
    print('break')
    print(spice.ktotal(KernelType.SPK.value))
    print(convert_UTC_to_ET(datetime.date.today()))
    frame, vector, number, bounds, bounds2 = find_fov(-226807, 4)

    print(bounds2)
