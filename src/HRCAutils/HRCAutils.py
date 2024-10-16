#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
from importlib_metadata import packages_distributions
from pathlib import Path

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-l', '--listets', is_flag=True, help='List entry points (i.e., utilities). (default: False)')
@click.version_option()
def main(listets):
	"""
Utilities of HRCA v1.0.

\b
Example:
  HRCAutils -l

\b
Date: 2024/10/16
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if listets:
		absdir=Path(__file__).parent
		scriptname=Path(__file__).stem
		for k,v in packages_distributions().items():
			if scriptname in v:
				print(','.join(v), k, sep='\t')

if __name__ == "__main__":
	main()
