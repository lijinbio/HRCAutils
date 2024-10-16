#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import shutil
import subprocess
import sys
import click
CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-v', '--verbose', is_flag=True, help='Verbose.')
@click.argument('cmdstr', nargs=-1)
def exeBash(cmdstr, condaenv, verbose):
	"""
Execute Bash commands.

\b
Example:
  exeBash date
  exeBash -e base date

\b
  exeBash which python3
  exeBash -e base which python3

\b
  # Calling original function
  from exeBash.exeBash import exeBash
  exeBash.callback(cmdstr=['which', 'python'], condaenv=None, verbose=False)
  exeBash.callback(cmdstr=['which', 'python'], condaenv='base', verbose=True)

\b
Note:
  1. Three parameters of exeBash.callback() are required. (No default values are allowed.)

  2. Revert back to stdin for a simplicity.

\b
  ```@2022/12/12
  2. Command strings are input from positional arguments. So, stdin is turned off by default for a security reason. (TODO)
  ```

\b
  ```
  -		if stdin:
  -			exitcode=os.system(f'bash -c "{cmdstr}"')
  -		else:
  -			exitcode=os.system(f'bash <<< "{cmdstr}"')
  +		exitcode=os.system(f'bash -c "{cmdstr}"')
  ```

\b
  3. Double quote are not allowed in the command string. Using the single quote instead. E.g.,

\b
  ```
  cmdstr="squeue -o '%.18i %.9P %.10j %.8u %.8T %.14M %.14l %.3C %.5D %.10m %R'"
  # cmdstr='squeue -o "%.18i %.9P %.10j %.8u %.8T %.14M %.14l %.3C %.5D %.10m %R"' # not this
  exeBash.callback(cmdstr=[cmdstr], condaenv=None, verbose=False)
  ```

\b
  4. Replace `os.system()` with `subprocess.call()`. @1/19/2023
  https://docs.python.org/3/library/subprocess.html?ref=hackernoon.com#replacing-os-system

\b
  5. Enable `bash` instead of /bin/sh by default for subprocess.call(shell=True). @2/10/2023

\b
  ```A test case for Bash feature
  exeBash 'date > >(md5sum -)'
  ```

\b
  6. The command prefix `set -euo pipefail;` is added to manage pipes. @3/28/2023

\b
  7. A micromamba environment is added to complement a conda environment for `-e|--condaenv`. @4/5/2023

\b
  8. A bug command using `conda/[micro]mamba` in `-e|--condaenv`. 
  The bash feature `set -o pipefail` will cause error in conda/[micro]mamba` commands. This is a bug. (TODO). @4/7/2023
  Thus, `set +euo pipefail` is necessary for `conda/[micro]mamba` commands.

\b
  ```E.g.,
  exeBash -e u_saturn 'micromamba list'; echo "$?" # will cause error
  exeBash -e u_saturn 'set +euo pipefail;micromamba list'; echo "$?" # `set +o pipefail` is necessary to turn off pipefail.
  ```

\b
Date: 2023/04/05
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	try:
		cmdprefix="set -euo pipefail;" # Adding prefix
		cmdstr=cmdprefix+' '.join(cmdstr)

		if condaenv:
			mamabrootprefix=os.environ.get('MAMBA_ROOT_PREFIX')
			if mamabrootprefix: # micromamba environment
				cmdstr=';'.join([
					f"source {mamabrootprefix}/etc/profile.d/micromamba.sh",
					f"micromamba activate {condaenv}",
					cmdstr,
					])
			else: # conda enviroment
				cmdstr=';'.join([
					'source "$(conda info --base)/etc/profile.d/conda.sh"',
					f"conda activate {condaenv}",
					cmdstr,
					])

		if verbose:
			click.echo(f'Start running: {cmdstr}', file=sys.stderr)

		# exitcode=os.system(f'bash -c "{cmdstr}"')
		# exitcode=subprocess.call(cmdstr, shell=True, executable=os.environ['SHELL'])
		exitcode=subprocess.call(cmdstr, shell=True, executable=shutil.which('bash'))

		if verbose:
			click.echo(f'Finish running: {cmdstr}', file=sys.stderr)

		if exitcode!=0:
			click.echo(f'Error: {cmdstr} failed. Exit code: {exitcode}', file=sys.stderr)
			sys.exit(exitcode)

	except OSError as e:
		click.echo("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)

	return exitcode

if __name__ == "__main__":
	sys.exit(exeBash())
