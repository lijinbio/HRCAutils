#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exePython.exePython import exePython
from pathlib import Path
import click
import random
import socket

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-k', '--batchkey', type=click.STRING, required=True, help='Batch key. E.g., sampleid.')
@click.option('-l', '--label', type=click.STRING, required=True, help='Cell type label.')
@click.option('-n', '--nhvg', type=click.INT, default=10000, show_default=True, help='Number of HVGs.')
@click.option('-f', '--flavor', type=click.Choice(['seurat', 'cell_ranger', 'seurat_v3']), is_flag=False, flag_value='seurat_v3', default='seurat', show_default=True, help='HVG algorithm.')
@click.option('-r', '--nlayer', type=click.INT, default=2, show_default=True, help='Number of hidden layers used for encoder and decoder NNs.')
@click.option('-t', '--nlatent', type=click.INT, default=30, show_default=True, help='Dimensionality of the latent space.')
@click.option('-p', '--epoch', type=click.INT, help='Max epoches for training.')
@click.option('-g', '--gpu', type=click.INT, default=-1, show_default=True, help='A GPU device.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, seed, batchkey, label, nhvg, flavor, nlayer, nlatent, epoch, gpu, infile):
	"""
To build a scArches reference using the scANVI model.

INFILE is a .h5ad file.

\b
Example:
  f=/storage/singlecell/jinli/wkfl/atlashumanprj/application/scArches/snRNA/BC/preproc/scrnah5adsubsetsamplingbykey/snRNA_BC.h5ad
  bname=$(basename "$f" .h5ad)
  outdir=$(mrrdir.sh)
  slurmtaco.sh -n g01 -- scarchesh5ad2refbyscanvi -d "$outdir" -b "$bname" -e scarches -k sampleid -l BC_macaque_name -n 10000 -- "$f"

\b
Note:
  1. A scANVI model might be always superior to a scVI model. So, use the scANVI model for now. (TODO: scVI model)
  2. This is forced to run on a GPU machine, see "-g|--gpu".
  3. HVG flavor.
  https://discourse.scverse.org/t/error-in-highly-variable-gene-selection/276/8

\b
  ```
  seurat v3 wants counts, the others want log normalized.
  ```

\b
See also:
  Depends:
    Python/scanpy
    Python/scvi-tools

\b
Date: 2023/08/22
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if gpu<0:
		if socket.gethostname()=='mhgcp-g00.grid.bcm.edu':
			gpu=random.randint(0, 1)
		elif socket.gethostname()=='mhgcp-g01.grid.bcm.edu':
			gpu=random.randint(0, 3)
		elif socket.gethostname() in ['gpu-8-0', 'gpu-8-1.mab', 'gpu-9-0']:
			gpu=0 # MAB GPU devices are managed by Slurm, so only device 0 is visable.
		else:
			click.echo(f"Error: please run on a GPU machine.", file=sys.stderr)
			sys.exit(-1)
	os.environ['CUDA_VISIBLE_DEVICES']=f"{gpu}"
	click.echo(f"Info: GPU device {os.environ['CUDA_VISIBLE_DEVICES']} is used.", file=sys.stderr)

	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"seed={seed}",
		f"batchkey='{batchkey}'",
		f"label='{label}'",
		f"nhvg={nhvg}",
		f"flavor='{flavor}'",
		f"nlayer={nlayer}",
		f"nlatent={nlatent}",
		f"epoch={epoch}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
