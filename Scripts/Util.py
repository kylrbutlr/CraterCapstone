"""
This module is for useful classes that may be implemented and segmented into their own modules.
"""

################
### lbl file ###
################

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


def get_lbl_value(lbl_text: str, var_name: str):
    """
    Returns the value of the varName that is declare in the lbl_text
    :param lbl_text: the contents of the lbl file
    :param var_name: the name of the variable declared in the lbl file
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


def get_image_time(lbl_text: str):
    """
    Returns the image time from the lbl file
    :param lbl_text: the contents of the lbl file
    :return: the time of the image, if the attribute isn't there return "0000-00-00"
    """

    time = get_lbl_value(lbl_text, "IMAGE_TIME")

    if time is None:
        time = "0000-00-00"

    return time


def get_image_size(lbl_text: str):
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