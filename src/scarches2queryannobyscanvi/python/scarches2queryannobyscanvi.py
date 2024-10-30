#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import scarches as sca
import seaborn as sns
from pathlib import Path
import anndata as ad

# reference model
labels_key='celltype'
batch_key='sampleid'

# query
if infile.endswith('.h5ad'):
	query=sc.read_h5ad(infile)
else:
	print(f'Error: format is not supported for {infile}')
	sys.exit(-1)

## query, labels_key
if labels_key in query.obs.columns:
	labels_key_original=f"{labels_key}_original"
	print(f'Warning: {labels_key=} exists in query, rename it to {labels_key_original}')
	query.obs.rename({labels_key: labels_key_original}, axis=1, inplace=True)
query.obs[labels_key]='Unknown'

## query, batch_key
if batchkey:
	if batchkey not in query.obs.columns:
		print(f'Error: {batchkey} is not in query. Double check the inputs.')
		sys.exit(-1)

	if batchkey!=batch_key: # copy a batch_key column for reference
		query.obs[batch_key]=query.obs[batchkey]
else:
	if batch_key not in query.obs.columns:
		print(f'Warning: {batch_key} used in reference is not in query. So, assume one single sample for query.')
		query.obs[batch_key]='_Query_scANVI_'

# prepare and load the query data
adata_query=sca.models.SCANVI.prepare_query_anndata(
	adata=query,
	reference_model=reference,
	inplace=False,
)
# surgery model
surgerymodel=sca.models.SCANVI.load_query_data(
	adata_query,
	reference,
	freeze_dropout=True,
)
print(f"Info: surgerymodel, labeled indices, {len(surgerymodel._labeled_indices)}")
print(f"Info: surgerymodel, unlabeled indices, {len(surgerymodel._unlabeled_indices)}")

# train
early_stopping_kwargs_surgery={
	'early_stopping_monitor': 'elbo_train',
	'early_stopping_patience': 10,
	'early_stopping_min_delta': 0.001,
	'plan_kwargs': {'weight_decay': 0.0},
	'check_val_every_n_epoch': 10,
}
surgerymodel.train(
	max_epochs=epoch,
	**early_stopping_kwargs_surgery
)
surgerymodel.save(f"{bname}_model_surgery", save_anndata=True)

# query latent
querylatent=sc.AnnData(
	surgerymodel.get_latent_representation(),
	obs=adata_query.obs,
)
querylatent.obs['scANVI_predictions']=surgerymodel.predict()
## predict probability
predprob=surgerymodel.predict(soft=True)
predprob['scANVI_prediction_max_probability']=predprob.max(axis=1)
predprob.insert(loc=0, column='barcode', value=predprob.index)
predprob.to_csv(f"{bname}_query_QC_predprob.txt.gz", sep='\t', index=False)
querylatent.obs['scANVI_prediction_max_probability']=predprob['scANVI_prediction_max_probability']

## visualization
sc.pp.neighbors(querylatent, random_state=seed)
sc.tl.umap(querylatent, random_state=seed)
sc.set_figure_params(dpi_save=500, figsize=(5, 5))
for splitby in querylatent.obs.columns:
	ncolor=len(querylatent.obs[splitby].value_counts())
	if ncolor<100:
		sc.pl.umap(querylatent, color=splitby, frameon=False, show=False, title='', save=f"{bname}_query_umap_{splitby}_wolabel.png")
		sc.pl.umap(querylatent, color=splitby, frameon=False, show=False, title='', save=f"{bname}_query_umap_{splitby}_ondata.png", legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal')
	else:
		palette=sns.husl_palette(ncolor)
		sc.pl.umap(querylatent, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_query_umap_{splitby}_wolabel.png")
		sc.pl.umap(querylatent, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_query_umap_{splitby}_ondata.png") # duplicate

## save the object
sc.write(filename=f'{bname}_query_latent.h5ad', adata=querylatent)
querylatent.obs.insert(loc=0, column='barcode', value=querylatent.obs.index)
querylatent.obs.to_csv(f'{bname}_query_latent_obs.txt.gz', sep='\t', index=False)
querylatent.var.insert(loc=0, column='symbol', value=querylatent.var.index)
querylatent.var.to_csv(f'{bname}_query_latent_var.txt.gz', sep='\t', index=False)
