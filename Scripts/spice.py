"""
This script declares the Spice function wrappers that are used in the project.
The functions in this script are renamed in order to improve readability.
"""

import spiceypy as spice
import enum
import numpy as np
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

###########
# Kernels #
###########


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


def get_shape_information(shape_file_path: str):
    """
    Loads the shape file (dsk) and returns the vertices of the shape and
    the vertices associated to each plate

    :param shape_file_path: relative path to the dsk files
    :return:
            Array of the coordinates of the vertices of the shape,
            Array of the indices of vertices that make up a plate
    """

    id = spice.dasopr(shape_file_path)
    handle = spice.dlabfs(id)

    num_vertices, num_plates = spice.dskz02(id, handle)

    vertices_array = spice.dskv02(id, spice.dlabfs(id), 1, num_vertices)
    plate_array = spice.dskp02(id, handle, 1, num_plates)

    return vertices_array, plate_array
