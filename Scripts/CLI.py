"""
This is the CLI module. More documentation to come.
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
    This a basic example a command line interface.
    This will be changed later on.

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
    
    click.echo("Spice Kernel File: %s\nExposure time: %s\nShape: %s\nAsteroid Image: %s" % (spice_kernel, ex_time, shape, ast_img))

if __name__ == '__main__':
    main()
