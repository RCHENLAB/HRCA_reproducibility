#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

# models
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import RadiusNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import MultinomialNB
from sklearn.neighbors import NeighborhoodComponentsAnalysis, KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.svm import SVC

# metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import scanpy as sc

# Model
if model=='GradientBoostingClassifier': # very slow
	clf=Pipeline([
		('scaler', StandardScaler()),
		('classifier', GradientBoostingClassifier(n_estimators=1024, random_state=seed, **params)),
		])
elif model=='KNeighborsClassifier':
	clf=Pipeline([
		('scaler', StandardScaler()),
		('classifier', KNeighborsClassifier(**params)),
		])
elif model=='NearestCentroid':
	clf=Pipeline([
		('scaler', StandardScaler()),
		('classifier', NearestCentroid(**params)),
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
		('classifier', RandomForestClassifier(n_estimators=1024, random_state=seed, **params)),
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

x=sc.read(infile)
X=x.X
y=x.obs[classlabel]
res=pd.DataFrame()
fold=0
skf=StratifiedKFold(n_splits=10)
for t_i, v_i in skf.split(X, y):
	fold+=1
	clf=clf.fit(X[t_i], y[t_i])
	pred=clf.predict(X[v_i])

	tab=pd.crosstab(y[v_i], pred)
	tab.to_csv(f'{bname}_tab_{fold}.txt.gz', sep='\t')

	tmp={
		'fold': fold,
		'accuracy': accuracy_score(y[v_i], pred),
		'precision': precision_score(y[v_i], pred, average='weighted'),
		'recall': recall_score(y[v_i], pred, average='weighted'),
		'f1': f1_score(y[v_i], pred, average='weighted'),
		# 'auc': roc_auc_score(y[v_i], clf.predict_proba(X[v_i]), average='weighted', multi_class='ovr') if model!='NearestCentroid' else -1,
		'balanced_acc': balanced_accuracy_score(y[v_i], pred),
		'ari': adjusted_rand_score(y[v_i], pred),
		}
	tmp=pd.DataFrame({k: [v] for k, v in tmp.items()})
	res=pd.concat([res, tmp], axis=0, join='outer', ignore_index=True)

res['fold']=res['fold'].apply(np.int64)
res.to_csv(f'{bname}_stats.txt.gz', sep='\t', index=False)
