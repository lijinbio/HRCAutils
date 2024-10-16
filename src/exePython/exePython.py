#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
from exeBash.exeBash import exeBash
import click
CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-e', '--expr', multiple=True, type=click.STRING, help='Expression. Multiple `-e|--expr` can be specified.')
@click.option('-s', '--script', type=click.Path(exists=True), help='An R script file to source.')
@click.option('-c', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-v', '--verbose', is_flag=True, help='Verbose.')
def exePython(expr, script, condaenv, verbose):
	"""
Execute Python expressions.

\b
Example:
  exePython -e "x=10; print(x)"
  exePython -e "x=10" -e "print(x)"
  exePython -e "x=10" -s script.py
  exePython -c base -v -e "print('expr: hello world.')" -e "x=10; print(x);" -s script.py
  ## exePython <<< 'print(10)' # stdin is not allowed

\b
  # calling original function
  from exePython.exePython import exePython
  exePython.callback(["x=10;print(x)"], script=None, condaenv=None, verbose=False)
  exePython.callback(["x=10", "print(x)"], script=None, condaenv='base', verbose=True)

\b
Note:
  1. Double quotation is needed for expressions (`-e|--expr`).

\b
See also:
  Related:
    exeR
  Depends:
    exeBash

\b
Date: 2022/12/05
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if script:
		expr+=(f"exec(open('{script}').read())",)
	expr=';'.join([e.rstrip(';') for e in expr])
	cmdstr=f'python3 -q <<< "{expr}"'
	return exeBash.callback([cmdstr], condaenv=condaenv, verbose=verbose)

if __name__ == "__main__":
	sys.exit(exePython())
