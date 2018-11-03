"""
This is the CLI module. More documentation to come.
"""
# Import packages/modules here
import click


@click.command()
@click.option('--count', default=1, help='Number of executions.')
@click.option('--spice_kernel', prompt='Spice kernel file: ', help='The file path of the spice kernel')
@click.option('--shape', default='Irregular', help='The configuration \
of the shape for the model')
def main(count, spice_kernel, shape):
    """
    This a basic example a command line interface.
    This will be changed later on.

    :param count: Number of executions.
    :type count: int
    :param spice_kernel: The file path of the spice kernel
    :type spice_kernel: str
    :param shape: The configuration of the shape for the model
    :type shape: str
    :return:
    """
    for x in range(count):
        click.echo('The shape is %s and the kernel is located at %s.' % (shape, spice_kernel))


if __name__ == '__main__':
    main()
