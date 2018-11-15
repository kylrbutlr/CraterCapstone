"""
This is the CLI module.

This is the input/output module for the system. It allows a user to
input required data (i.e. spice kernel file, asteroid image, shape, and execution time)
and generates a poly-line (list of coordinates) as output.
"""
# Import packages/modules here
import click

@click.command()
@click.option('--spice-kernel', prompt='Spice kernel file', help='The file path of the spice kernel')
@click.option('--ast-img', default='', help='Asteroid image')
@click.option('--shape', default='Irregular', help='The configuration of the shape for the model')
@click.option('--ex-time', default=0.0, help='Time of exposure')
def main(spice_kernel, shape, ast_img, ex_time):
    """
    This is the input/output of the system. To use it, simply run "python CLI.py" in a command window.
    Be sure that Python version 3.7.1 or higher is installed!

    CLI Version: 1.0

    :param spice-kernel: The file path of the spice kernel\n
    :type spice-kernel: str\n
    :param shape: The configuration of the shape for the model\n
    :type shape: str\n
    :param ast-img: Asteroid image\n
    :type shape: str\n
    :param ex-time: Time of exposure\n
    :type count: float\n
    :return: poly-line (array) of points that outline the visible portions of the specified body
    """

    output = [ (i, i + 1) for i in range(5) ] # TODO: replace with actual output after processing the input
    click.echo('Poly-line generated:')
    click.echo(output)

if __name__ == '__main__':
    main()
