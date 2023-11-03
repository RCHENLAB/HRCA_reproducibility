sam_hm=["5_S66", "6_S67", "7_S68", "8_S69"]
sam_mm=["1_S62", "2_S63", "3_S64", "4_S65"]

lane=["L001","L002","L003","L004"]

#/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/Enhancer_RNA_1_S1_L001_R1_001_mm_enhancer_count_reform_new

head="/storage/chenlab/Users/junwang/enhancer_validation/data/4-12-23/Enhancer_RNA_"
head1="/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/Enhancer_RNA_"

tail="R1_001_hm_enhancer_count_reform_new"

for sam in sam_hm:
	out=head + sam + "_" + tail
	ot=open(out, "w")
	count={}
	for la in lane:
		fn=head + sam + "_" + la +"_" + tail
		with open(fn, "r") as fn1:
			for line in fn1:
				info=line.rstrip().split()
				if info[0] not in count:
					count[info[0]]=0
				count[info[0]]+=int(info[1])
	sam_info=sam.split("_")
	sam_id=sam_info[0]+"_"+"S"+sam_info[0]
	fn_prev=head1 + sam_id + "_L001_R1_001_hm_enhancer_count_reform_new"
#	with open(fn_prev, "r") as fn1:
#		for line in fn1:
#			info=line.rstrip().split()
#			if info[0] not in count:
#				count[info[0]]=0
#			count[info[0]]+=int(info[1])
	for bc in count:
		ot.write(f'{bc}\t{count[bc]}\n')
	ot.close()


tail="R1_001_mm_enhancer_count_reform_new"

for sam in sam_mm:
	out=head + sam + "_" + tail
	ot=open(out, "w")
	count={}
	for la in lane:
		fn=head + sam + "_" + la +"_" + tail
		with open(fn, "r") as fn1:
			for line in fn1:
				info=line.rstrip().split()
				if info[0] not in count:
					count[info[0]]=0
				count[info[0]]+=int(info[1])
	sam_info=sam.split("_")
	sam_id=sam_info[0]+"_"+"S"+sam_info[0]
	fn_prev=head1 + sam_id + "_L001_R1_001_mm_enhancer_count_reform_new"
#	with open(fn_prev, "r") as fn1:
#		for line in fn1:
#			info=line.rstrip().split()
#			if info[0] not in count:
#				count[info[0]]=0
#			count[info[0]]+=int(info[1])
	for bc in count:
		ot.write(f'{bc}\t{count[bc]}\n')
	ot.close()

