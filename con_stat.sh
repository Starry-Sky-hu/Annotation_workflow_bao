#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J con_stat.sh
#SBATCH --qos=queue1
##################################################################
# @Author: huyong
# @Created Time : Tue Dec 14 10:30:34 2021

# @File Name: con_stat.sh
# @Description:
##################################################################

#for i in $(cat list)
#do
#	mkdir -p ${i}
#	cd ${i}
#	cp ../stat.sh .
#	sed -i 's/ref/'${i}'/g' stat.sh
#	cd ..
#done


# summary
for i in $(cat list)
do
	#echo -n "${i} " >> all_summary
	cat ${i}/summary/summary.txt >> all_summary
done

