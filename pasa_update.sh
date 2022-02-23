#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --qos=queue1
#SBATCH -w comput61
##################################################################
# @Author: huyong
# @Created Time : Tue Dec 14 20:50:44 2021

# @File Name: con.sh
# @Description:
##################################################################
threads=48

source /home/huyong/software/anaconda3/bin/activate mikado
cp /home/huyong/Yao_test/00_Ref/ref_genome.fa .
mkdir portcullis
cd portcullis
cp /home/huyong/Yao_test/ref_genome/02_Hisat/ref_genome_all.sort.bam .
portcullis prep -t ${threads} ../ref_genome.fa ref_genome_all.sort.bam
portcullis junc -t ${threads} ./portcullis_prep
portcullis filter -t ${threads} portcullis_prep ./portcullis_junc/portcullis.junctions.tab
cd ..

mkdir stringtie
cd stringtie
cp /home/huyong/Yao_test/ref_genome/02_Hisat/ref_genome_*2.sort.bam .
cd ..

source /home/huyong/software/anaconda3/bin/activate annotation_2

ls stringtie/*bam > bam_list
for i in $(cat bam_list)
do
	stringtie -p ${threads} -o ${i%.sort*}.gtf ${i}
done

source /home/huyong/software/anaconda3/bin/activate mikado

echo "#=====mikado_list.txt (only stringtie gtf)" > mikado_list.txt
ls stringtie/*gtf > gtf_list
cat gtf_list | sed 's/stringtie\///g' | sed 's/_hisat2.gtf//g' > sample_name
paste gtf_list sample_name | awk '{print $0"\tFalse\t1\tFalse\tTrue\tTrue"}' >> mikado_list.txt

mikado configure --list mikado_list.txt --reference ref_genome.fa --mode permissive --scoring plant.yaml --copy-scoring plant.yaml -bt ../uniprot_sprot_plants.fa --junctions portcullis/portcullis_filter/portcullis.pass.junctions.bed configuration.yaml

mikado prepare --json-conf configuration.yaml
diamond makedb --in ../uniprot_sprot_plants.fa --db uniprot -p ${threads}
diamond blastx --max-target-seqs 5 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop -q mikado_prepared.fasta -d uniprot -o mikado_prepared.blast.tsv -p ${threads}

singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit split2 -p ${threads} mikado_prepared.fasta -O mikado_prepared_fasta_split
cd mikado_prepared_fasta_split
ls *fasta > split_list

source /home/huyong/software/anaconda3/bin/activate annotation
for i in $(cat split_list)
do
	echo "#!/bin/bash" > ${i}.sh
	echo "#SBATCH --partition=big,low,smp01,smp02" >> ${i}.sh
	echo "#SBATCH -N 1" >> ${i}.sh
	echo "#SBATCH -c 1" >> ${i}.sh
	
	echo "/home/huyong/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ${i}" >> ${i}.sh

	echo "diamond blastp --max-target-seqs 1 --outfmt 6 --evalue 1e-5 -p 1 -q ${i}.transdecoder_dir/longest_orfs.pep -d ../uniprot > ${i}.transdecoder_dir/blastp.tsv" >> ${i}.sh

	echo "/home/huyong/software/interproscan/interproscan-5.53-87.0/bin/hmmer/hmmer3/3.1b1/hmmsearch --cpu 1 --domtblout ${i}.transdecoder_dir/pfam.domtblout /home/huyong/software/interproscan/interproscan-5.53-87.0/data/pfam/34.0/pfam_a.hmm ${i}.transdecoder_dir/longest_orfs.pep" >> ${i}.sh

	echo "/home/huyong/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t ${i} --retain_pfam_hits ${i}.transdecoder_dir/pfam.domtblout --retain_blastp_hits ${i}.transdecoder_dir/blastp.tsv" >> ${i}.sh
done

awk '{print "sh "$0".sh"}' split_list > split_parallel_jobs2run
parallel --jobs ${threads} < split_parallel_jobs2run
cd ..

cat mikado_prepared_fasta_split/mikado_prepared.part_*fasta.transdecoder.bed > mikado_prepared.fasta.transdecoder.bed
cat mikado_prepared_fasta_split/mikado_prepared.part_*fasta.transdecoder.cds > mikado_prepared.fasta.transdecoder.cds
cat mikado_prepared_fasta_split/mikado_prepared.part_*fasta.transdecoder.gff3 > mikado_prepared.fasta.transdecoder.gff3
cat mikado_prepared_fasta_split/mikado_prepared.part_*fasta.transdecoder.pep > mikado_prepared.fasta.transdecoder.pep

source /home/huyong/software/anaconda3/bin/activate mikado

mikado serialise -p ${threads} --json-conf configuration.yaml --xml mikado_prepared.blast.tsv --orfs mikado_prepared.fasta.transdecoder.gff3
mikado pick -p ${threads} --json-conf configuration.yaml --subloci-out mikado.subloci.gff3

singularity exec ~/contianer/OrthoGeneGL gffread -x mikado.loci.cds -y mikado.pep.fa -w mikado.cdna.fa -g ref_genome.fa mikado.loci.gff3

singularity exec /home/huyong/contianer/annotation_v0.1 \
	/home/software/PASApipeline-pasa-v2.5.1/bin/seqclean mikado.cdna.fa \
	-v /home/huyong/software/database/UniVec

singularity exec /home/huyong/contianer/annotation_v0.1 \
	/home/software/PASApipeline-pasa-v2.5.1/Launch_PASA_pipeline.pl \
	-c pasa.alignAssembly.config -C -R \
	-g ref_genome.fa \
	-t mikado.cdna.fa.clean -T \
	-u mikado.cdna.fa \
	--ALIGNERS blat,gmap \
	--CPU ${threads}

cp /home/huyong/Yao_test/ref_genome/05_maker/ref_genome.all.maker.filter0.75.gff3 .

singularity exec /home/huyong/contianer/annotation_v0.1 \
	/home/software/PASApipeline-pasa-v2.5.1/scripts/Load_Current_Gene_Annotations.dbi \
	-c pasa.alignAssembly.config -g ref_genome.fa -P ref_genome.all.maker.filter0.75.gff3

singularity exec /home/huyong/contianer/annotation_v0.1 \
	/home/software/PASApipeline-pasa-v2.5.1/Launch_PASA_pipeline.pl \
	-c pasa.annotationCompare.config -A \
	-g ref_genome.fa \
	-t mikado.cdna.fa.clean \
	--CPU ${threads}




#source /home/huyong/software/anaconda3/bin/activate annotation
#export PASAHOME=/home/huyong/software/PASApipeline-pasa-v2.5.1

#$PASAHOME/bin/seqclean mikado.cdna.fa -v /home/huyong/software/database/UniVec
#$PASAHOME/misc_utilities/accession_extractor.pl < mikado.cdna.fa.clean > pasa_fl_list

#$PASAHOME/Launch_PASA_pipeline.pl -c pasa.alignAssembly.config -C -R -g ref_genome.fa -f pasa_fl_list -t mikado.cdna.fa.clean -T -u mikado.cdna.fa --ALIGNERS gmap,minimap2

#cp /home/huyong/Yao_test/ref_genome/05_maker/ref_genome.all.maker.filter0.75.gff3 .

#$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi -c pasa.alignAssembly.config -g ref_genome.fa -P ref_genome.all.maker.filter0.75.gff3

#$PASAHOME/Launch_PASA_pipeline.pl -c pasa.annotationCompare.config -A -g ref_genome.fa -t mikado.cdna.fa.clean


