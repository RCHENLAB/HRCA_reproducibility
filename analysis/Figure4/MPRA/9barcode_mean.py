import pandas as pd
import sys

date=sys.argv[1]
spe=sys.argv[2]
eh_read={}

fn=f'/storage/chenlab/Users/junwang/enhancer_validation/data/{date}/enhancer_norm_{spe}-{date}'
#fn="/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/enhancer_norm_mm-3-24-23"

with open(fn,"r") as fn1:
	next(fn1)
	for line in fn1:
		info=line.rstrip().split()
		info1=info[0].split("_")
#		eh_id=info1[0]+"_"+info1[1]
		eh_id=info[0][0:-2:]
		print(eh_id)
		if info1[0] == "CRX":
			eh_id = "CRX"
		if eh_id not in eh_read:
			eh_read[eh_id]={}
		for i in range(4):
			if i not in eh_read[eh_id]:
				eh_read[eh_id][i]={}
				eh_read[eh_id][i]["num"]=int(0)
				eh_read[eh_id][i]["val"]=int(0)

			eh_read[eh_id][i]["num"] +=int(1)
			eh_read[eh_id][i]["val"] += float(info[i+1])

out=f'/storage/chenlab/Users/junwang/enhancer_validation/data/{date}/enhancer_norm_{spe}-{date}_enh'

#out="/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/enhancer_norm_mm-3-24-23_enh"


ot=open(out,"w")
#ot.write(f'\thm_1\thm_2\thm_3\thm_4\n')
ot.write(f'\trep_1\trep_2\trep_3\trep_4\n')

#ot.write(f'\tmm_1\tmm_2\tmm_3\tmm_4\n')
for eh_id in eh_read:
	ot.write(f'{eh_id}\t')
	for i in range(4):
		eh_read[eh_id][i]["mean"] = eh_read[eh_id][i]["val"] / eh_read[eh_id][i]["num"]
		ot.write(f'{eh_read[eh_id][i]["mean"]}\t')
	ot.write(f'\n')

ot.close() 	
		
