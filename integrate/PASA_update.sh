threads=48

source /home/huyong/software/anaconda3/bin/activate mikado
cp /home/huyong/grau_huyong/01_assembly/asm_results/Solanum_ochranthum_2hap.fa .
mkdir portcullis
cd portcullis
ln -s /home/huyong/grau_huyong/02_annotation_2/hisat2/Solanum_ochranthum.sort.bam* .
portcullis prep -t ${threads} ../Solanum_ochranthum_2hap.fa Solanum_ochranthum.sort.bam
portcullis junc -t ${threads} ./portcullis_prep
portcullis filter -t ${threads} portcullis_prep ./portcullis_junc/portcullis.junctions.tab
cd ..

mkdir stringtie
cd stringtie
ln -s /home/huyong/grau_huyong/02_annotation_2/stringtie/*gtf .
cd ..

source /home/huyong/software/anaconda3/bin/activate mikado

echo "#=====mikado_list.txt (only stringtie gtf)" > mikado_list.txt
ls stringtie/*gtf > gtf_list
cat gtf_list | sed 's/stringtie\///g' | sed 's/_hisat2.gtf//g' > sample_name
paste gtf_list sample_name | awk '{print $0"\tFalse\t1\tFalse\tTrue\tTrue"}' >> mikado_list.txt

mikado configure --list mikado_list.txt --reference Solanum_ochranthum_2hap.fa --mode permissive --scoring plant.yaml --copy-scoring plant.yaml -bt uniprot_sprot_plants.fa --junctions portcullis/portcullis_filter/portcullis.pass.junctions.bed configuration.yaml

mikado prepare --json-conf configuration.yaml
diamond makedb --in uniprot_sprot_plants.fa --db uniprot -p ${threads}
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

singularity exec ~/contianer/OrthoGeneGL gffread -x mikado.loci.cds -y mikado.pep.fa -w mikado.cdna.fa -g Solanum_ochranthum_2hap.fa mikado.loci.gff3

singularity exec /home/huyong/contianer/pasapipeline \
   /usr/local/src/PASApipeline/bin/seqclean mikado.cdna.fa \
   -v /home/huyong/software/database/UniVec_data/UniVec

singularity exec /home/huyong/contianer/pasapipeline \
   /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
   -c pasa.alignAssembly.config -C -R \
   -g Solanum_ochranthum_2hap.fa \
   -t mikado.cdna.fa.clean -T \
   -u mikado.cdna.fa \
   --ALIGNERS blat,gmap \
   --CPU ${threads}

cp /home/huyong/grau_huyong/02_annotation_2/maker/rename/filter_0.5/Solanum_ochranthum_2hap_filter_0.5.gff3 .

singularity exec /home/huyong/contianer/pasapipeline \
    /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
    -c pasa.alignAssembly.config -g Solanum_ochranthum_2hap.fa -P Solanum_ochranthum_2hap_filter_0.5.gff3

singularity exec /home/huyong/contianer/pasapipeline \
    /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
    -c pasa.annotationCompare.config -A \
    -g Solanum_ochranthum_2hap.fa \
    -t mikado.cdna.fa.clean \
    --CPU ${threads}

