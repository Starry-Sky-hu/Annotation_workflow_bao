#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -J con.sh
#SBATCH --qos=queue1
##################################################################
# @Author: huyong
# @Created Time : Tue Dec 14 20:50:44 2021

# @File Name: con.sh
# @Description:
##################################################################

for i in $(cat list)
do
	mkdir ${i}
	cd ${i}
	cp ../mikado.sh .
	cp ../pasa.alignAssembly.config .
	cp ../pasa.annotationCompare.config .
	sed -i 's/ref_genome/'${i}'/g' pasa.alignAssembly.config
	sed -i 's/ref_genome/'${i}'/g' pasa.annotationCompare.config
	sed -i 's/ref_genome/'${i}'/g' mikado.sh
	cd ..
done
