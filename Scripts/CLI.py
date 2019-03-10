"""
This is the CLI module.

This is the input/output module for the system. It allows a user to
input required data (i.e. spice kernel file, asteroid image, shape, and execution time)
and generates a poly-line (list of coordinates) as output.
"""
# Import packages/modules here
import click
import Geometry
"""
(body, nacid, utc, num_of_samples, target_body, target_body_reference, correction='NONE')
"""
@click.command()

# TODO add a separate option for each kernel file (Brian Vu and Nicole Darmawaskita)
@click.option('--spice-kernel', prompt='Spice kernel file', help='The file path of the spice kernel')
@click.option('--image', prompt='Image file', default='', help='Asteroid image')
# @click.option('--shape', default='Irregular', help='The configuration of the shape for the model')
@click.option('--samples', default=0, help='Number of samples to take from the image')
@click.option('--nacid', default=0, help='The nacid of the spacecraft')
@click.option('--target', default='', help='The target body that you want to observe')
@click.option('--reference', default='',
              help='The reference of the target body that you desire to observe')
@click.option('--spacecraft', default='', help='Instrument used to capture photo')
@click.option('--epoch', default=[], help='Time of exposure in utc')
@click.option('--filepath', default='../footprint.txt',
              help='The file path including the name of the output file')
def main(spice_kernel, spacecraft, image, epoch, samples, target, reference, nacid, filepath):
    """
    This is the input/output of the system. To use it, simply run "python CLI.py" in a command window.
    Be sure that Python version 3.7.1 or higher is installed!

    CLI Version: 1.0

    :param spice-kernel: The file path of the spice kernel\n
    :type spice-kernel: str\n
    :param shape: The configuration of the shape for the model\n
    :type shape: str\n
    :param image: Asteroid image\n
    :type shape: str\n
    :param instrument: Instrument used to capture photo\n
    :type instrument: str\n
    :param epoch: Time of exposure\n
    :type count: float\n
    :return: poly-line (array) of points that outline the visible portions of the specified body
    """
    Geometry.init_kernels()
    points = Geometry.Brute_Force(spacecraft, nacid, [epoch], samples, target, reference)
    output_Str = ""
    for point in points:
        output_Str += str(point[0]) + "," + str(point[1]) + "," + str(point[2]) + "\n"

    with open(filepath, 'w') as file:
        file.write(output_Str)

    """
    output = [ (i, i + 1) for i in range(5) ] # TODO: replace with actual output after processing the input
    click.echo('Poly-line generated:')
    click.echo(output)
    """


if __name__ == '__main__':
    main()
