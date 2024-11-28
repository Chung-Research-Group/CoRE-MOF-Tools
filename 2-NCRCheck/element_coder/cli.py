# -*- coding: utf-8 -*-

"""Command line interface for :mod:`element_coder`.

Why does this file exist, and why not put this in ``__main__``? You might be tempted to import things from ``__main__``
later, but that will cause problems--the code will get executed twice:

- When you run ``python3 -m element_coder`` python will execute``__main__.py`` as a script.
  That means there won't be any ``element_coder.__main__`` in ``sys.modules``.
- When you import __main__ it will get executed again (as a module) because
  there's no ``element_coder.__main__`` in ``sys.modules``.

.. seealso:: https://click.palletsprojects.com/en/7.x/setuptools/#setuptools-integration
"""


import click

from .data.coding_data import _PROPERTY_KEYS
from .decode import decode
from .encode import encode

__all__ = ("encode_cli", "decode_cli")


@click.command("cli")
@click.argument("element", type=str)
@click.option("--property", type=click.Choice(_PROPERTY_KEYS), default="mod_pettifor")
def encode_cli(element, property):
    """CLI for element_coder."""
    click.echo(encode(element, property))


@click.command("cli")
@click.argument("encoding", type=click.FLOAT, nargs=-1)
@click.option("--property", type=click.Choice(_PROPERTY_KEYS), default="mod_pettifor")
def decode_cli(encoding, property):
    """CLI for element_coder."""
    click.echo(decode(encoding, property))
