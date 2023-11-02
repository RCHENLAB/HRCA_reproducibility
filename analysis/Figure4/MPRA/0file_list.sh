#!/bin/sh
date=$1
ls /storage/chenlab/Users/junwang/enhancer_validation/data/${date}/*R1* | awk -F "." '{print $1}' | awk -F "/" '{print $NF}' > /storage/chenlab/Users/junwang/enhancer_validation/data/${date}/file_list


