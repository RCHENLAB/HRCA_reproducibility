#import scvelo as scv
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay, roc_curve, auc, roc_auc_score
from itertools import cycle
from sklearn.linear_model import LogisticRegression
from os.path import exists

indir=sys.argv[1]
cell=sys.argv[2]
label=sys.argv[3]

rna0=sc.read(f'{indir}/{cell}_rna-emb.h5ad')
atac0=sc.read(f'{indir}/{cell}_atac-emb.h5ad')

#dbt_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/emb_dbt.txt.gz"
dbt_table=f'{indir}/{cell}_emb_dbt.txt.gz'

rna=rna0
atac=atac0
rna_cell=rna0.obs.index
atac_cell=atac0.obs.index

rna0.obs.index=rna_cell
atac0.obs.index=atac_cell

if exists(dbt_table):
	dbt=pd.read_csv(dbt_table,sep="\t",index_col=0)
	rna_cell=rna0.obs.index.intersection(dbt.index)
	atac_cell=atac0.obs.index.intersection(dbt.index)
	rna=rna0[rna_cell,].copy()
	atac=atac0[atac_cell,].copy()


random_state=np.random.RandomState(1)

X=rna.obsm["X_glue"]
y=rna.obs[label]
X_atac=atac.obsm["X_glue"]

n_samples, n_features=X.shape
n_classes=len(np.unique(y))
#if n_classes==2 :
#	n_classes=1

target_names=np.unique(y)
##split the X_train and X_test into 50:50
(X_train, X_test, y_train, y_test) = train_test_split(X,y, test_size=0.5, stratify=y, random_state=1)

#######use SVM oneVsRestClassifier
#y_score_svm1=OneVsRestClassifier(LinearSVC(random_state=0)).fit(X_train, y_train).predict(X_test)
#y_score_svm=LabelBinarizer().fit_transform(y_score_svm1)

#y_score_svm_pred1=OneVsRestClassifier(LinearSVC(random_state=0)).fit(X, y).predict(X_atac)
#y_score_svm_pred=LabelBinarizer().fit_transform(y_score_svm_pred1)

#######use logistical regression
classifier= LogisticRegression(random_state=1, max_iter=1000)

#if n_classes==1:



y_score_lr=classifier.fit(X_train, y_train).predict_proba(X_test)
y_score_lr_pred=classifier.fit(X, y).predict_proba(X_atac)

#if n_classes==1 :
#	y_score_lr=y_score_lr[:,1]
#	y_score_lr_pred=y_score_lr_pred[:,1]
######binarize y and transform into one_hot_test
label_binarizer = LabelBinarizer().fit(y_train)
y_onehot_test = label_binarizer.transform(y_test)
y_onehot_test.shape 

####ROC curve using the OvR macro-average
# store the fpr, tpr, and roc_auc for all averaging strategies
fpr_lr, tpr_lr, roc_auc_lr = dict(), dict(), dict()
#fpr_svm, tpr_svm, roc_auc_svm = dict(), dict(), dict()
n=y_onehot_test.shape[1]
for i in range(n):
#	fpr_svm[i],tpr_svm[i], _ = roc_curve(y_onehot_test[:,i], y_score_svm[:,i])
#	roc_auc_svm[i] = auc(fpr_svm[i], tpr_svm[i])
	fpr_lr[i],tpr_lr[i], _ = roc_curve(y_onehot_test[:,i], y_score_lr[:,i])
	roc_auc_lr[i] = auc(fpr_lr[i], tpr_lr[i])



fpr_grid=np.linspace(0.0, 1.0, 1000)

#interpolate all ROC curves at these points
#mean_tpr_svm=np.zeros_like(fpr_grid)
mean_tpr_lr=np.zeros_like(fpr_grid)

for i in range(n_classes):
#	mean_tpr_svm += np.interp(fpr_grid, fpr_svm[i], tpr_svm[i])
	mean_tpr_lr += np.interp(fpr_grid, fpr_lr[i], tpr_lr[i])

#	mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])

# Average it an compute AUC 

#mean_tpr_svm /=n_classes

mean_tpr_lr /=n_classes


#fpr_svm["macro"] = fpr_grid
fpr_lr["macro"] = fpr_grid

#tpr_svm["macro"] = mean_tpr_svm
tpr_lr["macro"] = mean_tpr_lr

#roc_auc_svm["macro"] = auc(fpr_svm["macro"], tpr_svm["macro"])
roc_auc_lr["macro"] = auc(fpr_lr["macro"], tpr_lr["macro"])



fig, ax = plt.subplots(figsize=(6,6))

#plt.plot(
#	fpr_svm["macro"],
#	tpr_svm["macro"],
#	label=f'SVM macro-ave ROC curve (AUC = {roc_auc_svm["macro"]:.2f})',
#	color="navy",
#	linestyle=":",
#	linewidth=4,
#)


plt.plot(
	fpr_lr["macro"],
	tpr_lr["macro"],
	label=f'LR macro-ave ROC curve (AUC = {roc_auc_lr["macro"]:.2f})',
	color="red",
	linestyle=":",
	linewidth=4,
)


#colors =cycle(["aqua", "darkorange", "cornflowerblue"])
#for class_id, color in zip(range(n_classes), colors):
#	RocCurveDisplay.from_predictions(
#		y_onehot_test[:,class_id],
#		y_score_lr[:,class_id],
#		name=f'LR ROC curve for {target_names[class_id]}',
#		color=color,
#		ax=ax,
#	)

#colors=cycle(["red","blue","green"])
#colors=range(3)

#for class_id, color in zip(range(n_classes), colors):
#	RocCurveDisplay.from_predictions(
#		y_onehot_test[:,class_id],
#		y_score_svm[:,class_id],
#		name=f'SVM ROC curve for {target_names[class_id]}',
#		color=color,
#		ax=ax,
#	)

plt.plot([0,1],[0,1],"k--", label="ROC curve for chance level (AUC=0.5)")
plt.axis("square")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Extension of Receiver Operating Characteristic\nto One-vs-Rest multiclass")
plt.legend()
plt.savefig(f'{indir}/{cell}_auc.png')



#atac.obsm["svm_class"] =y_score_svm_pred
#df=pd.DataFrame(y_score_svm_pred, columns=range(y_score_svm_pred.shape[1]))
#atac.obs["svm_celltype"]=target_names[df.idxmax(axis=1)]

atac.obsm["lr_class"] =y_score_lr_pred
df=pd.DataFrame(y_score_lr_pred, columns=range(y_score_lr_pred.shape[1]))
atac.obs["lr_celltype"]=target_names[df.idxmax(axis=1)]

df1=df.max(axis=1)

df.to_csv(f'{indir}/{cell}_atac_obs_score.txt.gz',sep="\t")

atac.obs.to_csv(f'{indir}/{cell}_atac_obs.txt.gz',sep="\t")
atac.obs.iloc[df1[df1<0.9].index].to_csv(f'{indir}/{cell}_atac_obs_dbt.txt.gz',sep="\t")
#cm=atac.obs.groupby(["lr_celltype","svm_celltype"]).size().unstack(fill_value=0)
#cm.to_csv(f'{indir}/{cell}_lr_svm.txt.gz',sep="\t")

combined0=sc.read(f'{indir}/{cell}_combined-emb.h5ad')
#combined.obs["svm_celltype"]=[]
#combined.obs["lr_celltype"]=[]

new_cell=rna_cell.append(atac_cell)

combined=combined0[new_cell].copy()

#combined.obs.loc[atac.obs.index,"svm_celltype"]=atac.obs["svm_celltype"]
combined.obs.loc[atac.obs.index,"lr_celltype"]=atac.obs["lr_celltype"]

#combined.obs.loc[rna.obs.index,"svm_celltype"]=rna.obs[label]
combined.obs.loc[rna.obs.index,"lr_celltype"]=rna.obs[label]

combined=combined[combined.obs["lr_celltype"]!="NA"]

sc.pl.embedding(combined, basis="X_umap", color=["lr_celltype"],ncols=2,frameon=False,save=f'_{cell}_combined_pred.png') #, palette="tab20")

#sc.pl.embedding(combined, basis="X_umap", color=["svm_celltype","lr_celltype"],ncols=2,frameon=False,save=f'_{cell}_combined_pred.png') #, palette="tab20")
