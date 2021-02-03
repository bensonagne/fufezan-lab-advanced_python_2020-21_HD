#!/usr/bin/env python3
import sys
import click

@click.command()
@click.option('--name', default="???", help="Name to greet!")
def main(name="??"):
    """Print hello with input!

    Requires one parameter!!!
    """
    print(sys.argv)
    print("hello!", name)

if __name__ == '__main__':
    main()
    # if len(sys.argv) != 2:
    #     print(main.__doc__)
    # else:
    #     main(_input=sys.argv[1])
