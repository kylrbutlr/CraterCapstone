"""
This script executes the ray intercept functionality. Only instruments with rectangular field of views are implemented.
"""

from .spice import *
import numpy as np
from . import boundary
from spiceypy.utils.support_types import SpiceyError


def find_bound_extrema(fov_bounds):
    """
    Returns the min and max values of X and Y values of the fov bounds
    :param fov_bounds: the fov bounds of the instrument
    :return: minX: the minimum value of X values in fov bounds
                maxX: the maximum value of X values in fov bounds
                minY: the minimum value of Y values in fov bounds
                maxY: the maximum value of Y values in fov bounds
    """
    # Initialize values
    min_x = fov_bounds[0][0]
    max_x = fov_bounds[0][0]
    min_y = fov_bounds[0][1]
    max_y = fov_bounds[0][1]

    for i in range(len(fov_bounds)-1):
        x_value = fov_bounds[i + 1][0]
        y_value = fov_bounds[i + 1][1]

        if x_value > max_x:
            max_x = x_value
        if x_value < min_x:
            min_x = x_value

        if y_value > max_y:
            max_y = y_value
        if y_value < min_y:
            min_y = y_value

    return min_x, max_x, min_y, max_y


def get_direction_arrays(fov_bounds, num_of_samples):
    """
    Returns the direction vectors within the fovBounds
    :param fov_bounds: array of the edges of the field of view
    :param num_of_samples: the number of samples squared
    :return: direction_arrays: array of directional arrays
    """
    direction_arrays = []

    min_x, max_x, min_y, max_y = find_bound_extrema(fov_bounds)
    diff_x = max_x - min_x
    diff_y = max_y - min_y

    # Evenly distribute points across the fov bounds
    for i in range(num_of_samples):
        for j in range(num_of_samples):
            direction_arrays.append([min_x + i * (diff_x / num_of_samples), min_y + j * (diff_y / num_of_samples),
                                     fov_bounds[0][2]])

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
    reference_frame = fov[1]
    if fov[0] == 'RECTANGLE':
        direction_array = get_direction_arrays(fov[4], num_of_samples)

    for vec in direction_array:
        try:
            point, trgepc, srfvec = find_ray_surface_intercept_at_epoch(target_body, utc, target_body_reference,
                                                                        correction, body, reference_frame, vec)
            data.append(point)
        except SpiceyError:
            pass

    return data


def brute_force(body, nacid, utc, num_of_samples, target_body, target_body_reference, alpha_value):
    """
    Returns the footprint of the target body and the ray intercept points
    :param body: The observing or spacecraft body name
    :param nacid: The nacid value of the observing body
    :param utc: The epoch
    :param num_of_samples: The number of samples squared
    :param target_body: the naif name of the target body
    :param target_body_reference: the reference frame name of the target body
    :param alpha_value: the alpha value of the boundary
    :return: position_vectors: the coordinates of the edges of the footprint
                points: the collection of ray intercept points
    """

    position_vectors = calculate_position_vectors_of_ray_surface_intercept(body, nacid, utc, num_of_samples,
                                                                           target_body, target_body_reference)

    points = np.asarray(position_vectors)
    position_vectors = boundary.boundary(points, alpha_value)

    return position_vectors
