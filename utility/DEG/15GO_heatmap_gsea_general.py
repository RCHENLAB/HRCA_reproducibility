import os
import sys
cell=["Rod","Cone","BC","RGC","MG","AC","HC"]

mol=["KEGG","msigdbr","GO_BP"]

folder=sys.argv[1]
dem=sys.argv[2]
ct=int(sys.argv[3])
for m in mol:

	cell1={}
	cell1["up"]={}
	cell1["down"]={}
	go={}

	for c in cell:
#	fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/{c}_interval_cpm01_snRNA_rmH_FC1_young_cpm07_sim_rm0age_gender/{c}_GO_BP'
#		fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/{c}_interval_cpm01_snRNA_rmH_FC1_young_cpm07_sim_rm0age_gender/{c}_{m}'
		fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/{c}_interval_cpm01_snRNA_rmH_FC1_young_cpm07_{folder}/{c}_{m}'

#	fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/{c}_interval_cpm01_snRNA_rmH_FC1_young_cpm07_sim_rm0age_gender/{c}_KEGG'
		if os.path.exists(fn):
#		cell1[c]=1
			with open(fn,"r") as fn1:
				next(fn1)
				ct_down=0
				ct_up=0
				for line in fn1:
					info=line.split("\t")
					if float(info[8]) < 0.1:
						info[2]=info[2].replace("GOBP_","")
						label="up"
						if float(info[5]) <0:
							label="down"
							ct_down+=1
							cell1["down"][c]=1
						else:
							ct_up+=1
							cell1["up"][c]=1

						if label not in go:
							go[label]={}
				#print(info[2])
						if label=="up":
							if ct_up <=ct:
								if info[2] not in go[label]:
									go[label][info[2]]={}

								go[label][info[2]][c]=info[8]
						elif label=="down":
							if ct_down <=ct:

								if info[2] not in go[label]:
									go[label][info[2]]={}


								go[label][info[2]][c]=info[8]


	for l in go:
#	fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top15_gender_{l}'  # msigdbr top 10
		fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top{ct}_{m}_{dem}_{l}'  # msigdbr top 10

#		fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top10_{m}_gender_{l}'  # msigdbr top 10
#	fn=f'/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top15_KEGG_gender_{l}'  # KEGG top 10

		with open(fn, "w") as fn1:
			line=[]
			for c in cell1[l]:
				line.append(str(c))
		
			line1="\t".join(line)
			fn1.write(f'{line1}\n')

			for key in go[l]:
				line=str(key)
				for c in cell1[l]:
					if c in go[l][key]:
						line+="\t"+go[l][key][c]

					else:
						line+="\t"+"1"

				fn1.write(f'{line}\n')



						
