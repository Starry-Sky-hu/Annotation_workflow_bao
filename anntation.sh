
# fastp 
R1=($(cat raw_R1_list))
R2=($(cat raw_R2_list))
for i in `seq 0 2`
do
	fastp -a auto \
		--adapter_sequence_r2 auto \
		--detect_adapter_for_pe \
		-w 10 \
		-i ${R1[i]} -I ${R2[i]} \
		-o ${R1[i]%.fq.gz*}_clean.fq.gz -O ${R2[i]%.fq.gz*}_clean.fq.gz \
		--n_base_limit 5 \
		--cut_window_size 4 \
		--cut_mean_quality 20 \
		--length_required 75 \
		--qualified_quality_phred 15
done


# mapping
hisat2-build -p 52 C151_final_H1_chr_EDTA_rm.fa ./C151_final_H1_chr_EDTA_rm

hisat2 --dta \
	-x ../ref/C151_final_H1_chr_EDTA_rm \
	-1 ../rawdata/C151_R1.clean.fq.gz \
	-2 ../rawdata/C151_R2.clean.fq.gz \
	--rna-strandness RF \
	--summary-file {output.summary} \
	--new-summary \
	-p 52 \
	-S C151_final_H1_chr_EDTA_rm_hisat2.sam

samtools view -@ 52 -Sb C151_final_H1_chr_EDTA_rm_hisat2.sam > C151_final_H1_chr_EDTA_rm_hisat2.bam
samtools sort -@ 52 -o C151_final_H1_chr_EDTA_rm_hisat2_sorted.bam C151_final_H1_chr_EDTA_rm_hisat2.bam

# braker
set +u;source /home/huyong/software/anaconda3/bin/activate annotation; set -u

braker.pl \
	--species=C151_H1_braker \
	--genome=../ref/C151_final_H1_chr_EDTA_sm.fa \
	--bam=../hisat_result/C151_final_H1_chr_EDTA_rm_hisat2_sorted.bam \
	--workingdir=. \
	--gff3 \
	--nocleanup \
	--softmasking \
	--cores 48


cd ../maker
/home/huyong/software/mpich/bin/mpiexec -n 52 maker -base C151_H1 maker_opts.ctl maker_exe.ctl maker_bopts.ctl

# stringtie

set +u;source /home/huyong/software/anaconda3/bin/activate annotation_2; set -u

stringtie -p 52 --rf -o C151_final_H1.gtf ../hisat_result/C151_final_H1_chr_EDTA_rm_hisat2_sorted.bam

gffread -E C151_final_H1.gtf -o - | sed "s#transcript#match#g" | sed "s#exon#match_part#g" > C151_final_H1.gff3

# maker
source ~/software/anaconda3/bin/activate annotation

/home/huyong/software/mpich/bin/mpiexec -n 52 maker -base C151_H1 maker_opts.ctl maker_exe.ctl maker_bopts.ctl
gff3_merge -d C151_H1.maker.output/C151_H1_master_datastore_index.log
fasta_merge -d C151_H1.maker.output/C151_H1_master_datastore_index.log

# rename
for i in $(cat list)
do
	mkdir -p ${i}
	cd ${i}
	echo "cp /home/huyong/Yao_test/${i}/05_maker/${i}.all.maker.gff3 ." > rename.sh
	echo "perl ~/software/maker_AED_rename_S.pl -n ${i}_ -gff ${i}.all.maker.gff3" >> rename.sh
	echo "mv AED_1.maker.gff ${i}.all.rename.gff3" >> rename.sh
	echo "perl ~/software/maker/bin/quality_filter.pl -a 0.75 ${i}.all.rename.gff3 > ${i}.all.rename_0.75.gff3" >> rename.sh
	echo "perl /home/wuyaoyao/03-Solanaceae/src/longestCDS_PrimaryGff.pl ${i}.all.rename_0.75.gff3 > ${i}.all.rename_0.75_best.gff3" >> rename.sh
	echo "/home/wuyaoyao/software/iTools_Code/bin/iTools Gfftools getCdsPep -Ref ~/Yao_test/00_Ref/${i}.fa -Gff ${i}.all.rename_0.75_best.gff3 -OutPut ${i}" >> rename.sh
	echo "gunzip *gz" >> rename.sh
	cd ..
done


for i in $(cat list)
do
	cd ${i}
	sh rename.sh
	sleep 5
	cd ..
done

