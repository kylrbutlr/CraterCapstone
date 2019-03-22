"""
This is the CLI module.

This is the input/output module for the system. It allows a user to
input required data (i.e. spice kernel file, asteroid image, shape, and execution time)
and generates a poly-line (list of coordinates) as output.
"""
# Import packages/modules here
import click
import Geometry
import pandas as pd
import Spice
import Util
"""
(body, nacid, utc, num_of_samples, target_body, target_body_reference, correction='NONE')
"""
@click.command()

# TODO add a separate option for each kernel file (Brian Vu and Nicole Darmawaskita)
@click.option('--kernel_path', prompt='Spice meta-kernel file', help='The file path of the spice meta-kernel')
@click.option('--lbl_path', prompt='Image lbl file', default='', help='Asteroid image')
# @click.option('--shape', default='Irregular', help='The configuration of the shape for the model')
@click.option('--observing_body', prompt='The observing body or spacecraft instrument', default='',
              help='Instrument used to capture photo')
@click.option('--target_body', prompt='Target Body', default='', help='The target body that you want to observe')
@click.option('--target_body_frame',prompt='Target Body Reference Frame',  default='',
              help='The reference of the target body that you desire to observe')

def main(kernel_path:str, lbl_path:str, observing_body:str, target_body:str, target_body_frame:str):
    """
    This is the input/output of the system. To use it, simply run "python CLI.py" in a command window.
    Be sure that Python version 3.7.1 or higher is installed!

    CLI Version: 1.0

    :param kernelPath: the path of the metakernel
    :param lblPath: the path of the lbl file
    :param observingBody: the naif name of the observing body
    :param targetBody: the naif name of the target body
    :return: an array of the end points of the edges of the footprint
    """
    Spice.load_kernel(kernel_path)
    time, size = Util.get_lbl_information(lbl_path)
    size = int(size ** (1 / 2))
    footprint, points = Geometry.brute_force(observing_body, Spice.get_id_code(observing_body),
                                             Spice.convert_utc_to_et(time), size, target_body, target_body_frame)

    footprint.shape = (footprint.shape[0], footprint.shape[1], footprint.shape[2])
    pd.DataFrame(points).to_csv('../sampleFootprint.csv', header=None, index=None)

    Spice.unload_all_kernels()

    """
    output = [ (i, i + 1) for i in range(5) ] # TODO: replace with actual output after processing the input
    click.echo('Poly-line generated:')
    click.echo(output)
    """

if __name__ == '__main__':
    main()

