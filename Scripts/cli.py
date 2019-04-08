"""
This is the Command Line Interface module.

This is the input/output module for the system. It allows a user to
input required data and generates a poly-line (list of coordinates) as output.
"""

# Import packages/modules here
import click
import geometry
import pandas as pd
import spice
import util


@click.command()
@click.option('--kernel_path', prompt='Spice meta-kernel file', help='The file path of the spice meta-kernel')
@click.option('--lbl_path', prompt='Image lbl file', default='', help='Asteroid image')
# @click.option('--shape', default='Irregular', help='The configuration of the shape for the model')
@click.option('--observing_body', prompt='The observing body or spacecraft instrument', default='',
              help='Instrument used to capture photo')
@click.option('--target_body', prompt='Target Body', default='', help='The target body that you want to observe')
@click.option('--target_body_frame', prompt='Target Body Reference Frame',  default='',
              help='The reference of the target body that you desire to observe')
@click.option('--sample_size_scalar', prompt='Sample Size Scalar',  default='5',
              help='The integer scalar of samples that is multiplied by the size of the image file')
@click.option('--boundary_alpha_value', prompt='Boundary Alpha Value',  default='.15',
              help='The boundary alpha value that determines the accuracy of the footprint')
@click.option('--file_location', prompt='Output File Path',  default='',
              help='The path of the output file')
@click.option('--file_name', prompt='Output File name',  default='output',
              help='The output file name')
def main(file_name: str, file_location:str, kernel_path: str, lbl_path: str, observing_body: str, target_body: str,
         target_body_frame: str, sample_size_scalar: int, boundary_alpha_value: float):
    """
    This is the input/output of the system. To use it, simply run "python CLI.py" in a command window.
    Be sure that Python version 3.7.1 or higher is installed!

    CLI Version: 1.0

    :param file_name: the name of the output file
    :param file_location: the path location of the output file
    :param kernel_path: the path of the metakernel
    :param lbl_path: the path of the lbl file
    :param observing_body: the naif name of the observing body
    :param target_body: the naif name of the target body
    :param target_body_frame: the reference frame of the target body
    :param sample_size_scalar: the scalar of the amount of samples taken
    :param boundary_alpha_value: the alpha value for determining the footprint boundary
    :return: an array of the end points of the edges of the footprint
    """

    spice.load_kernel(kernel_path)
    time, size = util.get_lbl_information(lbl_path)
    size = int((size ** (1 / 2)) * int(sample_size_scalar))

    footprint = geometry.brute_force(observing_body, spice.get_id_code(observing_body), spice.convert_utc_to_et(time),
                                     size, target_body, target_body_frame, float(boundary_alpha_value))

    polyline = util.convert_to_polyline(footprint)
    formatted_polyline = util.prepare_to_save_to_file(polyline)
    pd.DataFrame(formatted_polyline).to_csv(file_location + file_name + '.csv', header=None, index=None)

    spice.unload_all_kernels()


if __name__ == '__main__':
    main()

