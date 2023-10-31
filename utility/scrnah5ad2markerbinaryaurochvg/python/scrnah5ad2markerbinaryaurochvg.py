#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

# models
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import RadiusNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import MultinomialNB
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier # XGBClassifier

# metrics
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import RocCurveDisplay
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# main
import matplotlib.pyplot as plt
from multiprocess import Pool
import itertools
import sys
import pandas as pd
import scanpy as sc
import jlutilspy as jl

def get_classifier(model, seed, params):
	if model=='GradientBoostingClassifier': # very slow
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', GradientBoostingClassifier(random_state=seed, **params)),
			])
	elif model=='XGBClassifier': # a fast version of XGBoost
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', XGBClassifier(random_state=seed, **params)),
			])
	elif model=='KNeighborsClassifier':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', KNeighborsClassifier(**params)),
			])
	elif model=='RadiusNeighborsClassifier':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', RadiusNeighborsClassifier(**params)),
			])
	elif model=='LogisticRegression':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', LogisticRegression(multi_class='ovr', random_state=seed, **params))
			])
	elif model=='MultinomialNB':
		clf=Pipeline([
			('classifier', MultinomialNB(**params)),
			])
	elif model=='NeighborhoodComponentsAnalysis': # very slow
		clf=Pipeline([
			('scaler', StandardScaler()),
			('nca', NeighborhoodComponentsAnalysis(random_state=seed, **params)),
			('knn', KNeighborsClassifier(**params)),
			])
	elif model=='RandomForestClassifier':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', RandomForestClassifier(random_state=seed, **params)),
			])
	elif model=='SGDClassifier':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', SGDClassifier(loss='log', random_state=seed, **params)),
			])
	elif model=='SVC':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', SVC(kernel='rbf', random_state=seed, probability=True, decision_function_shape='ovr', **params)),
			])
	elif model=='SVClinear':
		clf=Pipeline([
			('scaler', StandardScaler()),
			('classifier', SVC(kernel='linear', random_state=seed, probability=True, decision_function_shape='ovr', **params)),
			])
	return clf

def get_classifier_noscale(model, seed, params):
	if model=='GradientBoostingClassifier': # very slow
		clf=Pipeline([
			('classifier', GradientBoostingClassifier(random_state=seed, **params)),
			])
	elif model=='XGBClassifier': # a fast version of XGBoost
		clf=Pipeline([
			('classifier', XGBClassifier(random_state=seed, **params)),
			])
	elif model=='KNeighborsClassifier':
		clf=Pipeline([
			('classifier', KNeighborsClassifier(**params)),
			])
	elif model=='RadiusNeighborsClassifier':
		clf=Pipeline([
			('classifier', RadiusNeighborsClassifier(**params)),
			])
	elif model=='LogisticRegression':
		clf=Pipeline([
			('classifier', LogisticRegression(multi_class='ovr', random_state=seed, **params))
			])
	elif model=='MultinomialNB':
		clf=Pipeline([
			('classifier', MultinomialNB(**params)),
			])
	elif model=='NeighborhoodComponentsAnalysis': # very slow
		clf=Pipeline([
			('nca', NeighborhoodComponentsAnalysis(random_state=seed, **params)),
			('knn', KNeighborsClassifier(**params)),
			])
	elif model=='RandomForestClassifier':
		clf=Pipeline([
			('classifier', RandomForestClassifier(random_state=seed, **params)),
			])
	elif model=='SGDClassifier':
		clf=Pipeline([
			('classifier', SGDClassifier(loss='log', random_state=seed, **params)),
			])
	elif model=='SVC':
		clf=Pipeline([
			('classifier', SVC(kernel='rbf', random_state=seed, probability=True, decision_function_shape='ovr', **params)),
			])
	elif model=='SVClinear':
		clf=Pipeline([
			('classifier', SVC(kernel='linear', random_state=seed, probability=True, decision_function_shape='ovr', **params)),
			])
	return clf

def list2combn(lst, n):
	return list(itertools.combinations(lst, n))

def cmd(*args):
	print(f"{args=}", file=sys.stderr)
	group, feature=args
	c_x=x[:, feature].copy()

	if verbose:
		print(f"{c_x=}, {feature=}")

	X=c_x.X
	y=c_x.obs[label]
	y=(y==group).astype(int).tolist()

	if verbose:
		print(f"{group=}, {X=}")
		print("y:")
		print(jl.npflat2valuecounts(y))

	del c_x # delete to save memory

	c_clf=clf.fit(X, y)
	pred=c_clf.predict(X)

	if verbose:
		print("pred:")
		print(jl.npflat2valuecounts(pred))
		print("crosstab: y, pred")
		tab=pd.crosstab(y, pred)
		tab.to_csv(sys.stdout, sep='\t')
		precision=precision_score(y, pred, average='binary')
		print(f"{precision=}")

	tmp={
		'group': group,
		'symbol': ','.join(feature),
		'precision': precision_score(y, pred, average='binary'),
		'recall': recall_score(y, pred, average='binary'),
		'f1': f1_score(y, pred, average='binary'),
		'auc': roc_auc_score(y, c_clf.predict_proba(X)[:, 1]), # binary, see https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html
		}
	tmp=pd.DataFrame({k: [v] for k, v in tmp.items()})

	if verbose:
		tmp.to_csv(sys.stderr, sep='\t')

		imgfile=f"{outdir}/{bname}_{group}_{'.'.join(feature)}.pdf"
		plt.figure(figsize=(4, 4), dpi=500)
		RocCurveDisplay.from_predictions(y, pred)
		plt.tight_layout()
		plt.savefig(imgfile)
		plt.close()

	return tmp

def scanpyobj2topgenes(adata, group, nhvg, ngene, outfile):
	sc.pp.normalize_total(adata)
	sc.pp.log1p(adata)
	sc.pp.highly_variable_genes(
		adata,
		flavor='seurat',
		n_top_genes=nhvg,
		subset=True,
		batch_key=batchkey,
	)
	# sc.tl.rank_genes_groups(adata, groupby=group, n_genes=ngene, method='wilcoxon')
	sc.tl.rank_genes_groups(adata, groupby=group, n_genes=None, method='wilcoxon')

	groups=adata.uns['rank_genes_groups']['names'].dtype.names
	header=['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
	result=pd.DataFrame()
	for g in groups:
		tmp=pd.concat(
			[
				pd.DataFrame(adata.uns['rank_genes_groups'][h])[g]
				for h in header
			],
			axis=1,
			keys=header,
			)
		tmp.insert(loc=0, column='group', value=g)
		result=pd.concat([result, tmp], axis=0, join='outer', ignore_index=True)
	result=result[(result['logfoldchanges']>=log2foldchange) & (result['pvals_adj']<=qvalue) & (result['scores']>0)]
	result.sort_values(by=['group', 'scores'], axis=0, ascending=[True, False], inplace=True, ignore_index=True, key=None)
	result=result.groupby('group').head(ngene)
	result.to_csv(outfile, sep='\t', index=False)
	topgenes=dict(result.groupby('group')['names'].apply(list))
	return topgenes

if __name__ == '__main__':
	x=sc.read(infile)
	sc.pp.filter_genes(x, min_counts=1) # raw counts, sum of counts per gene
	print(f"{x=}")
	print(x.obs[label].value_counts())
	
	if topgenefile:
		topgenes=pd.read_csv(topgenefile, sep='\t', header=0)
		topgenes=dict(topgenes.groupby('group')['names'].apply(list))
	else:
		topgenes=scanpyobj2topgenes(x.copy(), group=label, nhvg=nhvg, ngene=ngene, outfile=f"{outdir}/{bname}_topgenes.txt.gz")
	print(f"Info: {topgenes=}")

	if group:
		topgenes={k: topgenes[k] for k in group}
		print(f"Info: {group=}; {topgenes=}")

	if norm:
		sc.pp.normalize_total(x)
		sc.pp.log1p(x)

	clf=get_classifier_noscale(model, seed, params)
	
	result=pd.DataFrame()
	for group, symbols in topgenes.items():
		# run the job in parallel
		with Pool(processes=numthreads) as p:
			greslist=p.starmap(cmd, itertools.product([group], list2combn(symbols, ncombn)))
		result=pd.concat([result]+greslist)
		print(f"Info: {group=}, {result=}")

	result.sort_values(by=['group', 'auc'], axis=0, ascending=[True, False], inplace=True, ignore_index=True, key=None)
	result.to_csv(f"{outdir}/{bname}.txt.gz", sep='\t', index=False)
