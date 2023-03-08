threads=52
mkdir AED_filter_0.75 AED_filter_0.5 all

cp /home/huyong/grau_huyong/02_annotation_2/PASA_update/filter_0.75/Solanum_ochranthum_hap1_PASA_db_filter_0.75.gene_structures_post_PASA_updates.gff3 AED_filter_0.75
cp /home/huyong/grau_huyong/02_annotation_2/PASA_update/filter_0.5/Solanum_ochranthum_hap1_PASA_db_filter_0.5.gene_structures_post_PASA_updates.gff3 AED_filter_0.5
cp /home/huyong/grau_huyong/02_annotation_2/PASA_update/all/Solanum_ochranthum_hap1_gene_structures_post_PASA_updates.gff3 all 

ln -s ../../../01_assembly/asm_results/Solanum_ochranthum_hap1.fa .

ls */*gff3 > gff3_list
for i in $(cat gff3_list)
do
	singularity exec ~/contianer/OrthoGeneGL \
		gffread ${i} \
		-g Solanum_ochranthum_hap1.fa \
		-y ${i%.gff3*}_protein.fa \
		-w ${i%.gff3*}_exon.fa \
		-x ${i%.gff3*}_cds.fa \
		-S
done

#############   mapping ratio  ##################
ls */*exon.fa > exon_list
ln -s ../../fastp .
ls fastp/*/*gz > fq_list
fq=($(cat fq_list))

source activate annotation_2
for exon in $(cat exon_list)
do
	kallisto index -i ${exon}.idx ${exon}
done

let max=${#fq[@]}-1
for exon in $(cat exon_list)
do
	for i in `seq 0 2 ${max}`
	do
		t=${fq[i]%clean*}
		out=${t##*/}
		kallisto quant -i ${exon}.idx -o ${exon%/*}/${out}kallisto_result -t ${threads} ${fq[i]} ${fq[i+1]}
	done
done

#############   busco     #################
ls */*protein.fa > protein_list
for i in $(cat protein_list)
do
   cd ${i%/*}
   singularity exec /home/huyong/contianer/busco_container \
       busco -m protein \
           -i ${i#*/} \
           -o ${i#*/}_busco \
           -l /home/huyong/software/solanales_odb10 \
           -c ${threads} \
           --offline
   cd ..
done

#############    pfam   #################
source activate java11

for i in $(cat protein_list)
do
	/home/huyong/software/interproscan/interproscan-5.55-88.0/interproscan.sh -appl PfamA \
	   -iprlookup -goterms -f tsv \
	   -dp -cpu ${threads} \
	   -b ${i%.fa*}.pfam \
	   -i ${i}
done


source activate R_Python
for i in $(cat gff3_list)
do
	cd ${i%/*}
	python -m jcvi.annotation.stats stats ${i#*/} &> jcvi_stats 
	python -m jcvi.annotation.stats genestats ${i#*/} &> jcvi_genestats
	python -m jcvi.annotation.stats histogram ${i#*/}
	python -m jcvi.annotation.stats statstable ${i#*/} &> jcvi_statstable
	python -m jcvi.annotation.stats summary ${i#*/} ../Solanum_ochranthum_hap1.fa &> jcvi_summary
 	cd ..
done


mkdir -p summary
#echo -n "all pfam: " > summary.txt
cat all/*.pfam.tsv | cut -f1 | sort | uniq | wc -l > summary/pfam_summary.txt
#echo -n "AED_filter_0.75 pfam: " >> summary.txt
cat AED_filter_0.75/*.pfam.tsv | cut -f1 | sort | uniq | wc -l >> summary/pfam_summary.txt
#echo -n "AED_filter_0.5 pfam: " >> summary.txt
cat AED_filter_0.5/*.pfam.tsv | cut -f1 | sort | uniq | wc -l >> summary/pfam_summary.txt

cat all/*_busco/short* | grep "C:" | sed 's/\t//g' | sed 's/ //g' > summary/busco_summary.txt
cat AED_filter_0.75/*_busco/short* | grep "C:" | sed 's/\t//g' | sed 's/ //g' >> summary/busco_summary.txt
cat AED_filter_0.5/*_busco/short* | grep "C:" | sed 's/\t//g' | sed 's/ //g' >> summary/busco_summary.txt

grep "Gene" all/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' > summary/gene_summary.txt
grep "Gene" AED_filter_0.75/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/gene_summary.txt
grep "Gene" AED_filter_0.5/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/gene_summary.txt

grep "Exon" all/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' > summary/exon_summary.txt
grep "Exon" AED_filter_0.75/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/exon_summary.txt
grep "Exon" AED_filter_0.5/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/exon_summary.txt

grep "Intron" all/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' > summary/intron_summary.txt
grep "Intron" AED_filter_0.75/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/intron_summary.txt
grep "Intron" AED_filter_0.5/jcvi_summary | awk '{print $6,$4,$5,$7,$3}' >> summary/intron_summary.txt

grep "p_pseudoaligned" all/*kallisto_result/run_info.json | awk '{print $3}' | sed 's/,//g' | xargs -n100 > summary/mapping_ratio
grep "p_pseudoaligned" AED_filter_0.75/*kallisto_result/run_info.json | awk '{print $3}' | sed 's/,//g' | xargs -n100 >> summary/mapping_ratio
grep "p_pseudoaligned" AED_filter_0.5/*kallisto_result/run_info.json | awk '{print $3}' | sed 's/,//g' | xargs -n100 >> summary/mapping_ratio


paste summary/busco_summary.txt summary/pfam_summary.txt summary/mapping_ratio summary/gene_summary.txt summary/exon_summary.txt summary/intron_summary.txt > summary/summary.txt

sed -i 's/ /_/g' summary/summary.txt




