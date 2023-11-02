import sys

enhancer=sys.argv[1] #"enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_hm_bc"
date=sys.argv[4]
######map_enhancer="/storage/chenlab/Users/junwang/enhancer_validation/data/1-22-23/Hu_enhancer_JX_S81_L003_R1_001_hm_enhancer_count_reform_new"
#map_enhancer="/storage/chenlab/Users/junwang/enhancer_validation/data/2-5-23/JX7_Enhancer_hg_CRX_hm_enhancer_count_reform_new"
#map_enhancer=f'/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/{sys.argv[2]}_{sys.argv[3]}_enhancer_count_reform_new'
map_enhancer=f'/storage/chenlab/Users/junwang/enhancer_validation/data/{date}/{sys.argv[2]}_{sys.argv[3]}_enhancer_count_reform_new'

############enhancer="enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_mm_bc"
#map_enhancer="/storage/chenlab/Users/junwang/enhancer_validation/data/1-22-23/Ms_enhancer_JX_S82_L003_R1_001_mm_enhancer_count_reform_new"
###############map_enhancer="/storage/chenlab/Users/junwang/enhancer_validation/data/2-5-23/JX8_Enhancer_ms_mm_enhancer_count_reform_new"

#>IRD_chr1:2049486-2049720
eh_orig={}
eh_seq={}
with open(enhancer, "r") as eh:
	for line in eh:
		if line[0]==">":
			ID=line.strip().replace(">","")
			info=line.strip().replace(">","").split("_") 
			if info[0] not in eh_orig:
				eh_orig[info[0]]={}

			eh_orig[info[0]][ID]=1
			eh_seq[ID]=next(eh)
			#next(eh)


eh_map={}
with open(map_enhancer, "r" ) as eh:
	for line in eh:

#QTL_chr8:455058-457469_QTL      17      0.000191756717125003    88654   3.86948456189468e-05    439335
		ID=line.strip().replace(">","").split("\t")
		if int(ID[1]) > 0:
			info=line.strip().replace(">","").split("_")
			if info[0] not in eh_map:
				eh_map[info[0]]={}
			eh_map[info[0]][ID[0]]=1

ot=f'{map_enhancer}_sum'

out=open(ot,"w")

for tp in eh_orig:
	count=len(eh_orig[tp])
	count1=0
	if tp in eh_map:
		count1=len(eh_map[tp])

	out.write(f'{tp}\t{count}\t{count1}\n')

ot1=f'{map_enhancer}_sum_drop'

out1=open(ot1,"w")
	
for tp in eh_orig:
	for bc in eh_orig[tp]:
		if bc not in eh_map[tp]:
			out1.write(f'{tp}\t{bc}\n')


