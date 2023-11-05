import sys

label_sign={}
sep="&"
#fn=["Currant_36848389PG_IST","Currant_36848389PG_ONL","Currant_36848389PG_OST","FritscheLG2016_26691988NG","Gharahkhani_33627673NC","Hysi_32231278NG","Jiang_34737426NG"]
fn=["all_TRUE_hg19"]
for fn1 in fn:
#	file1=f'/storage/chenlab/Users/junwang/human_meta/data/finemap/{fn1}_pip_var_all_merged_uniform_anno'
	file1=f'/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/{fn1}_retina_gencodeHg19_1000g_flt_pip_var_merged_new_anno'
#	file1=f'/storage/chenlab/Users/junwang/human_meta/data/age_DESeq2_batch_age_dream/{ct}_interval_cpm01_snRNA_rmH/{ct}_DEG_res'
	with open (file1,"r") as f1:
		for line in f1:
			info=line.rstrip().split()
			if info[-3] != "anno;":
				#anno;3_UTR|peak|
#				feature=info[-3].replace("anno;","").replace("3_UTR","Exon").replace("5_UTR","Exon").split("|")
				feature=info[-3].replace("anno;","").split("|")

				ft_line=sep.join(feature)
				if ft_line not in label_sign:
					label_sign[ft_line]=0
				label_sign[ft_line]+=1



#output=f'/storage/chenlab/Users/junwang/human_meta/data/finemap/pip_var_all_merged_uniform_anno_all'
output=f'/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/all_TRUE_hg19_retina_gencodeHg19_1000g_flt_pip_var_merged_new_anno_sum'
out=open(output,"w") 

arr=[]
for label in label_sign:
	num=label_sign[label]
#	out.write(f'"{label}"={num},')
	arr.append(f'"{label}"={num}')

new_line="\n".join(arr)

out.write(f'{new_line}')

out.close()

