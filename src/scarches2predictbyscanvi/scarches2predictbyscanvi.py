#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exePython.exePython import exePython
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-r', '--reference', type=click.Path(exists=True, resolve_path=True), required=True, help='A pre-trained model directory for reference.')
@click.option('-k', '--batchkey', type=click.STRING, help='Batch key in query. E.g., sampleid.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-p', '--epoch', type=click.INT, help='Max epoches for training.')
@click.option('-g', '--gpu', type=click.INT, default=0, show_default=True, help='A GPU device.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, reference, batchkey, seed, epoch, gpu, infile):
	"""
To predict cell type label against a pre-trained scArches reference using the scANVI model.

INFILE is a query .h5ad file.

\b
Example:
  reference=HRCA_snRNA_allcells_model_scanvi
  infile=Retina_sample1.h5ad
  bname=$(basename "$infile" .h5ad)
  outdir=$(mrrdir.sh)
  scarches2predictbyscanvi -d "$outdir" -b "$bname" -e scarches -r "$reference" -k sampleid -- "$infile"

\b
Note:
  1. The batch key (by "-k|--batchkey") is for query, and it might be different from the batch_key used in the reference. Internally, the reference batch_key will be added to the query.

\b
See also:
  Upstream:
    scarchesh5ad2refbyscanvi
  Depends:
    Python/scarches
    Python/scanpy
    Python/scvi-tools

\b
Date: 2024/10/16
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if gpu>0:
		os.environ['CUDA_VISIBLE_DEVICES']=f"{gpu}"
	click.echo(f"Info: GPU device {os.environ['CUDA_VISIBLE_DEVICES']} is used.", file=sys.stderr)

	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"reference='{reference}'",
		f"batchkey='{batchkey}'" if batchkey is not None else f"batchkey=None",
		f"seed={seed}",
		f"epoch={epoch}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
