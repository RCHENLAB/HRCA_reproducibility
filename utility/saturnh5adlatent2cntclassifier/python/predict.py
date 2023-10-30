#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, adjusted_rand_score
import scanpy as sc
import jlutilspy as jl
import pandas as pd

x=sc.read(infile)
X=x.X
model=jl.cpickleload(model)
pred=model.predict(X)

res=pd.DataFrame({
	'barcode': x.obs.index,
	'predict': pred,
	})
res.to_csv(f'{bname}.txt.gz', sep='\t', index=False)

if label:
	y=x.obs[label]
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
