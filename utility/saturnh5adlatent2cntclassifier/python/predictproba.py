#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, adjusted_rand_score
import scanpy as sc
import jlutilspy as jl
import pandas as pd

x=sc.read(infile)
X=x.X
model=jl.cpickleload(model)
predproba=model.predict_proba(X)

res=pd.DataFrame(predproba, columns=model.classes_)
res['prob_max']=predproba.max(axis=1)
res['predict_max']=model.classes_[predproba.argmax(axis=1)]
res['predict']=res['predict_max']
res.loc[res['prob_max']<threshold, 'predict']='Unassigned'
res.insert(loc=0, column='barcode', value=x.obs.index)
res.to_csv(f'{bname}.txt.gz', sep='\t', index=False)

if label:
	y=x.obs[label]
	tabmax=pd.crosstab(y, res['predict_max'].tolist())
	tabmax.to_csv(f'{bname}_max.txt.gz', sep='\t')

	pred=res['predict'].tolist()
	tab=pd.crosstab(y, pred)
	tab.to_csv(f'{bname}_tab.txt.gz', sep='\t')

	tmp={
		'accuracy': accuracy_score(y, pred),
		'precision': precision_score(y, pred, average='weighted'),
		'recall': recall_score(y, pred, average='weighted'),
		'f1': f1_score(y, pred, average='weighted'),
		## Not used: labels of y are not necessarily the same as pred.
		# 'auc': roc_auc_score(y, model.predict_proba(X), average='weighted', multi_class='ovr'),
		'balanced_acc': balanced_accuracy_score(y, pred),
		'ari': adjusted_rand_score(y, pred),
		}
	res=pd.DataFrame({k: [v] for k, v in tmp.items()})
	res.to_csv(f'{bname}_stats.txt.gz', sep='\t', index=False)
