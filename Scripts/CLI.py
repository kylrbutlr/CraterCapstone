import click

"""
This is an example of using click for a command line interface. 
"""
"""@click.command()
@click.option('--count', default=1, help='Number of greetings.')
@click.option('--name', prompt='Your name', help='The person to greet.')
def hello(count, name):
    \"""Simple program that greets NAME for a total of COUNT times.\"""
    for x in range(count):
        click.echo('Hello %s!' % name)
"""


@click.command()
@click.option('--count', default=1, help='Number of greetings.')
@click.option('--spice_kernel', prompt='Spice kernel file: ', help='The file path of the spice kernel')
@click.option('--Shape', default='Irregular', prompt='Configuration of shape: ', help='The configuration \
of the shape for the model')
def main(count, spice_kernel, shape):
    """
    This is the main method for the command line interface.
    Although the main method takes in no parameters, the options
    from the command line will act as parameters.

    Will be refactored and modified later

    :param: count, spice_kernel, shape
    :return: None

    """
    for x in range(count):
        click.echo('The shape is %s and the kernel is located at %s.' % (shape, spice_kernel))


if __name__ == '__main__':
    main()
