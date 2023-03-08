## TE annotation
source activate EDTA
perl /home/huyong/software/EDTA-master/EDTA.pl --genome Solanum_ochranthum_2hap.fa --anno 1 --sensitive 1 -t 52

bedtools maskfasta -fi Solanum_ochranthum_2hap.fa -bed Solanum_ochranthum_2hap.fa.mod.EDTA.TEanno.gff3 -fo Solanum_ochranthum_2hap_hardmask.fa
bedtools maskfasta -soft -fi Solanum_ochranthum_2hap.fa -bed Solanum_ochranthum_2hap.fa.mod.EDTA.TEanno.gff3 -fo Solanum_ochranthum_2hap_softmask.fa



## fastp clean RNAseq reads
fastp --detect_adapter_for_pe \
  -w 2 \
  -i ~/grau_huyong/SeqData/RNAseq/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_1.fq.gz \
  -I ~/grau_huyong/SeqData/RNAseq/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_2.fq.gz \
  -o ~/grau_huyong/02_annotation/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_clean_1.fq.gz \
  -O ~/grau_huyong/02_annotation/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_clean_2.fq.gz \
  --json ~/grau_huyong/02_annotation/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_fastp.json \
  --html ~/grau_huyong/02_annotation/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_fastp.html



## hisat2 mapping to hardmask genome
hisat2-build -p 52 Solanum_ochranthum_2hap_hardmask.fa Solanum_ochranthum_2hap_hardmask

source activate annotation_2
hisat2 --dta -x Solanum_ochranthum_2hap_hardmask \
  -1 /home/huyong/grau_huyong/02_annotation_2/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_clean_1.fq.gz \
  -2 /home/huyong/grau_huyong/02_annotation_2/fastp/Solanum_ochranthum_ML_1/Solanum_ochranthum_ML_1_clean_2.fq.gz \
  --summary-file Solanum_ochranthum_ML_1_clean_hisat_summary \
  --new-summary -p 13 > Solanum_ochranthum_ML_1_clean.sam
samtools sort -@ 13 -O bam -o Solanum_ochranthum_ML_1_clean_sort.bam Solanum_ochranthum_ML_1_clean.sam

## merge multiple samples
bam=`ls *sort.bam`
echo ${bam}
samtools merge -@ 52 -O bam -o Solanum_ochranthum.bam ${bam}
samtools sort -@ 52 -o Solanum_ochranthum.sort.bam Solanum_ochranthum.bam
samtools index Solanum_ochranthum.sort.bam



## stringtie convert BAM format to gff3 format
stringtie -p 13 -o Solanum_ochranthum_ML_1_clean.gtf Solanum_ochranthum_ML_1_clean_sort.bam
gffread -E Solanum_ochranthum_ML_1_clean.gtf -o - | sed "s#transcript#match#g" | sed "s#exon#match_part#g" > Solanum_ochranthum_ML_1_clean.gff3



## braker
## NOTE: cores not permit more than 48 !!!
source activate annotation
braker.pl \
        --species=Solanum_ochranthum_phase_edta_sensitive \
        --genome=Solanum_ochranthum_2hap_softmask.fa \
        --bam=Solanum_ochranthum.sort.bam \
        --workingdir=. \
        --gff3 \
        --nocleanup \
        --softmasking \
        --cores 48



## maker
source activate annotation
mpiexec -n 52 maker -base Solanum_ochranthum_2hap

fasta_merge -d Solanum_ochranthum_2hap.maker.output/Solanum_ochranthum_2hap_master_datastore_index.log
gff3_merge -n -g -o Solanum_ochranthum_2hap_ng.gff3 -d Solanum_ochranthum_2hap.maker.output/Solanum_ochranthum_2hap_master_datastore_index.log


