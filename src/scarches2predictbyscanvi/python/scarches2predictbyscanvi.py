#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import scarches as sca
import seaborn as sns
from pathlib import Path
import anndata as ad

# reference model
refmodel=sca.models.SCANVI.load(
	reference,
	# use_gpu=True, # removed from scVI 1.1.0
	)
labels_key=refmodel.registry_['setup_args']['labels_key']
batch_key=refmodel.registry_['setup_args']['batch_key']
print(f"Info: refmodel, labeled indices, {len(refmodel._labeled_indices)}")
print(f"Info: refmodel, unlabeled indices, {len(refmodel._unlabeled_indices)}")

# query
if infile.endswith('.h5ad'):
	query=sc.read(infile)
else:
	print(f'Error: format is not supported for {infile}')
	sys.exit(-1)

## query, labels_key
if labels_key in query.obs.columns:
	labels_key_original=f"{labels_key}_original"
	print(f'Warning: {labels_key=} exists in query, rename it to {labels_key_original}')
	query.obs.rename({labels_key: labels_key_original}, axis=1, inplace=True)
query.obs[labels_key]=refmodel.unlabeled_category_

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
	# use_gpu=True, # removed from scVI 1.1.0
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
	# use_gpu=True, # deleted from scVI 1.1.0
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

# co-embedding between ref and query
refadatafile=f"{reference}/adata.h5ad"
queryadatafile=f"{bname}_model_surgery/adata.h5ad"
if Path(refadatafile).exists() and Path(queryadatafile).exists():
	refadata=sc.read(refadatafile)
	queryadata=sc.read(queryadatafile)
	fulldata=ad.concat(adatas={'reference': refadata, 'query': queryadata}, axis=0, join='outer', label='scANVI_dataset', merge=None)
	# co-embedding
	coemb=sc.AnnData(
		surgerymodel.get_latent_representation(adata=fulldata),
		obs=fulldata.obs,
	)
	coemb.obs['scANVI_predictions']=surgerymodel.predict(adata=fulldata)
	
	## predict probability
	predprob=surgerymodel.predict(adata=fulldata, soft=True)
	predprob['scANVI_prediction_max_probability']=predprob.max(axis=1)
	predprob.insert(loc=0, column='barcode', value=predprob.index)
	predprob.to_csv(f"{bname}_coemb_QC_predprob.txt.gz", sep='\t', index=False)
	coemb.obs['scANVI_prediction_max_probability']=predprob['scANVI_prediction_max_probability']
	
	## visualization
	sc.pp.neighbors(coemb, random_state=seed)
	sc.tl.umap(coemb, random_state=seed)
	sc.set_figure_params(dpi_save=500, figsize=(5, 5))
	for splitby in coemb.obs.columns:
		ncolor=len(coemb.obs[splitby].value_counts())
		if ncolor<100:
			sc.pl.umap(coemb, color=splitby, frameon=False, show=False, title='', save=f"{bname}_coemb_umap_{splitby}_wolabel.png")
			sc.pl.umap(coemb, color=splitby, frameon=False, show=False, title='', save=f"{bname}_coemb_umap_{splitby}_ondata.png", legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal')
		else:
			palette=sns.husl_palette(ncolor)
			sc.pl.umap(coemb, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_coemb_umap_{splitby}_wolabel.png")
			sc.pl.umap(coemb, color=splitby, palette=palette, frameon=False, show=False, title='', save=f"{bname}_coemb_umap_{splitby}_ondata.png") # duplicate
	
	## save the object
	sc.write(filename=f'{bname}_coemb_latent.h5ad', adata=coemb)
	coemb.obs.insert(loc=0, column='barcode', value=coemb.obs.index)
	coemb.obs.to_csv(f'{bname}_coemb_latent_obs.txt.gz', sep='\t', index=False)
	coemb.var.insert(loc=0, column='symbol', value=coemb.var.index)
	coemb.var.to_csv(f'{bname}_coemb_latent_var.txt.gz', sep='\t', index=False)

