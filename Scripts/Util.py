"""
This module is for useful classes that may be implemented and segmented into their own modules.
"""

################
### lbl file ###
################

import spiceypy as sp
import numpy as np
def get_lbl_information(path: str):
    """
    Returns the image time and size from the lbl file
    :param path: path of the lbl file
    :return: image time formatted as a date as a string
                image size in bytes
    """

    lbl_text = load_image_lbl(path)
    time = get_image_time(lbl_text)
    size = get_image_size(lbl_text)

    return time, size

def load_image_lbl(path: str):
    """
    Loads the lbl file and returns the text
    :param path: path of the lbl file
    :return: text in the lbl file
    """

    lbl_text = open(path, "r").read()
    return lbl_text

def get_lbl_value(lbl_text:str, var_name: str):
    """
    Returns the value of the varName that is declare in the lbl_text
    :param lbl_text: the contents of the lbl file
    :param varName: the name of the variable declared in the lbl file
    :return: the value of the variable name in as a string
    """

    var_name = var_name.upper()
    text_iterator = iter(lbl_text.splitlines())
    while True:
        try:
            line = next(text_iterator)

            if var_name in line:
                return line.split("= ")[1].replace(" ", "")

        except StopIteration:
            break

def get_image_time(lbl_text:str):
    """
    Returns the image time from the lbl file
    :param lbl_text: the contents of the lbl file
    :return: the time of the image, if the attribute isn't there return "0000-00-00"
    """

    time = get_lbl_value(lbl_text, "IMAGE_TIME")

    if time is None:
        time = "0000-00-00"

    return time

def get_image_size(lbl_text:str):
    """
    Returns the size of the image from the lbl file
    :param lbl_text: the contents of the lbl file
    :return: the size of the image, if the attribute isn't there then return 0
    """

    size = get_lbl_value(lbl_text, "FILE_RECORDS")
    if size is not None:
        size = int(size)
    else:
        size = 0

    return size


def convert_to_polyline(points):
    point_dict = dict()
    for point1, point2 in points:
        if tuple(point1) in point_dict:
            continue
        else:
            point_dict[tuple(point1)] = tuple(point2)
    result = []
    while point_dict:
        key = next(iter(point_dict.values()))
        val = point_dict[key]
        polyline_points = [point_dict.pop(key)]
        while val != key:
            polyline_points.append(np.asarray(val))
            val = point_dict.pop(val)
        polyline_points.append(np.asarray(val))  # Remove this line if line does not need first point

        result.append(polyline_points)
    result = np.asarray(result)
    return result


def convert_points_to_lat_long(points):
    lat_long_points = []
    for point in points:
        radius, lat, long = sp.reclat(point)
        # lat = lat*sp.dpr() # converts to plantocentric coordinates
        # long = long*sp.dpr() # converts to plantocentric coordinates
        lat_long_points.append([lat, long, radius])
    lat_long_points = np.asarray(lat_long_points)
    return lat_long_points


def correct_output(points):
    polyline_result = []
    polyline = convert_to_polyline(points)
    for line in polyline:
        polyline_result.append(convert_points_to_lat_long(line))
    return polyline_result


def prepare_to_save_to_file(polyline):
    formatted_polylines = []
    for poly in polyline:
        formatted_polyline = []
        for point in poly:
            points_string = str(point[0]) + " " + str(point[1]) + " " + str(point[2])
            formatted_polyline.append(points_string)
        formatted_polylines.append(formatted_polyline)
    formatted_polylines = np.asarray(formatted_polylines)
    return formatted_polylines

