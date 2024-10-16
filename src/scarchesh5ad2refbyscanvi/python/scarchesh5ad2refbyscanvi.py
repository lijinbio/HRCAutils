#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import scarches as sca
import seaborn as sns
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, adjusted_rand_score

if infile.endswith('.h5ad'):
	x=sc.read_h5ad(infile)
else:
	print(f'Error: format is not supported for {infile}')
	sys.exit(-1)

# bug fix filtering genes with zero count
gene_subset, _=sc.pp.filter_genes(x, min_counts=1, inplace=False)
x=x[:, gene_subset].copy()

# raw counts
x.layers['hvgcounts']=x.X.copy() # for raw counts and will be subset by HVG
# subset by HVGs
## seurat_v3 use raw counts, and others use normalized counts
if flavor!='seurat_v3':
	sc.pp.normalize_total(x)
	sc.pp.log1p(x)
sc.pp.highly_variable_genes(
	x,
	flavor=flavor,
	n_top_genes=nhvg,
	subset=True,
	batch_key=batchkey,
	)

# raw count matrix
xdata=sc.AnnData(
	x.layers['hvgcounts'],
	obs=x.obs,
	var=x.var,
	)
del x ## free x to save memory

# scVI
sca.models.SCVI.setup_anndata(
	xdata,
	batch_key=batchkey,
	labels_key=label,
	)
# see scvi.module.MULTIVAE
# https://docs.scvi-tools.org/en/stable/api/reference/scvi.module.MULTIVAE.html
vae=sca.models.SCVI(
	xdata,
	n_layers=nlayer,
	n_latent=nlatent,
	use_batch_norm='none',
	use_layer_norm='both',
	deeply_inject_covariates=False,
	encode_covariates=False,
)
vae.train(
	max_epochs=epoch,
	# use_gpu=True, # deleted in scVI v1.1.0
	)
vae.save(f"{bname}_model_scvi", save_anndata=True)

# scANVI
sca.models.SCANVI.setup_anndata(
	xdata,
	batch_key=batchkey,
	labels_key=label,
	unlabeled_category='Unknown',
	)
scanvae=sca.models.SCANVI.from_scvi_model(
	vae,
	unlabeled_category='Unknown',
	)
print(f"Info: Labelled Indices: {len(scanvae._labeled_indices)}")
print(f"Info: Unlabelled Indices: {len(scanvae._unlabeled_indices)}")

scanvae.train(
	max_epochs=epoch,
	# use_gpu=True, # deleted in scVI v1.1.0
	)
scanvae.save(f"{bname}_model_scanvi", save_anndata=True)

xlatent=sc.AnnData(
	scanvae.get_latent_representation(),
	obs=xdata.obs,
	)
xlatent.obs['scANVI_predictions']=scanvae.predict()

## training statistics
y=xlatent.obs[label]
ypred=xlatent.obs['scANVI_predictions']
tab=pd.crosstab(y, ypred)
tab.to_csv(f"{bname}_QC_tab.txt.gz", sep='\t')
res=pd.DataFrame([{
	'accuracy': accuracy_score(y, ypred),
	'precision': precision_score(y, ypred, average='weighted'),
	'recall': recall_score(y, ypred, average='weighted'),
	'f1': f1_score(y, ypred, average='weighted'),
	'balanced_acc': balanced_accuracy_score(y, ypred),
	'ari': adjusted_rand_score(y, ypred),
	}])
res.to_csv(f"{bname}_QC_stats.txt.gz", sep='\t', index=False)
## prediction probabilities
predprob=scanvae.predict(soft=True)
predprob.insert(loc=0, column='barcode', value=predprob.index)
predprob.to_csv(f"{bname}_QC_predprob.txt.gz", sep='\t', index=False)

## visualization
sc.pp.neighbors(xlatent, random_state=seed)
sc.tl.umap(xlatent, random_state=seed)
sc.set_figure_params(dpi_save=500, figsize=(5, 5))
for splitby in xlatent.obs.columns:
	ncolor=len(xlatent.obs[splitby].value_counts())
	if ncolor<100:
		sc.pl.umap(xlatent, color=splitby, frameon=False, show=False, title='', save=f"{bname}_umap_{splitby}_wolabel.png")
		sc.pl.umap(xlatent, color=splitby, frameon=False, show=False, title='', save=f"{bname}_umap_{splitby}_ondata.png", legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal')
	else:
		palette=sns.husl_palette(ncolor)
		sc.pl.umap(xlatent, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_umap_{splitby}_wolabel.png")
		sc.pl.umap(xlatent, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_umap_{splitby}_ondata.png") # duplicate

# save the object
sc.write(filename=f'{bname}_latent.h5ad', adata=xlatent)
xlatent.obs.insert(loc=0, column='barcode', value=xlatent.obs.index)
xlatent.obs.to_csv(f'{bname}_obs.txt.gz', sep='\t', index=False)
xlatent.var.insert(loc=0, column='symbol', value=xlatent.var.index)
xlatent.var.to_csv(f'{bname}_var.txt.gz', sep='\t', index=False)

