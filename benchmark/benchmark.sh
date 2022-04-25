#########################
# Section One
#########################

# Download the test dataset.
wget -c https://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_fastqs.tar

# Unpack it
tar xvf neurons_900_fastqs.tar

# Download cell barcode candiate list
wget -c https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

# Convert FASTQ to FASTQ+
/usr/bin/time -v PISA parse2 -rule 'CR,R1:1-16,3M-february-2018.txt.gz,CB,1;UR,R1:17-26;R1,R2' -1 reads.fq -report parse.csv neurons_900_fastqs/neurons_900_S1_L001_R1_001.fastq.gz,neurons_900_fastqs/neurons_900_S1_L002_R1_001.fastq.gz neurons_900_fastqs/neurons_900_S1_L001_R2_001.fastq.gz,neurons_900_fastqs/neurons_900_S1_L002_R2_001.fastq.gz -nw

# Perform alignment
STAR --genomeDir star --readFilesIn reads.fq --outSAMunmapped Within --runThreadN 20 --outStd SAM 1>aln.sam

# Convert SAM to BAM and adjust mapq for multihits
/usr/bin/time -v PISA sam2bam -adjust-mapq -gtf genes.gtf aln.sam -t 4 -@ 4 -o aln.bam -report aln.csv 2>sam2bam.log

# Gene annotation
/usr/bin/time -v PISA anno -gtf genes.gtf aln.bam -o anno.bam -t 4 -report anno.csv 2>anno.log


# UMI correction
/usr/bin/time -v PISA corr -tag UR -new-tag UB -tags-block CB,GN -o corr.bam anno.bam -@ 4 2>corr.log

# Sort BAM
sambamba sort corr.bam -o sorted.bam -t 4

# Remove PCA duplicates for each molecular, consider Cell barcodes and UMIs
/usr/bin/time -v PISA rmdup -tags CB,UR -@ 4 -o rmdup.bam sorted.bam -nw 2>rmdup.log

# Variant calling
bcftools mpileup -Ou -f ~/Documents/datasets/GRCh38/fasta/genome.fasta ./rmdup.bam  | bcftools call -vmO z -o var.vcf.gz

# Annotate genetic variants
/usr/bin/time -v PISA anno -vtag VF -vcf var.vcf.gz rmdup.bam -o anno_vcf.bam -t 4 2>anno_vcf.log

# Count genotype X cell matrix
mkdir VF
/usr/bin/time -v PISA count -cb CB -anno-tag VF -umi UB -outdir var -@ 4 anno_vcf.bam 2>count_var.log    


#########################
# Section Two
#########################

wget -c https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fastqs.tar

tar xvf atac_pbmc_1k_nextgem_fastqs.tar

# Parse read 1 and read 2 
/usr/bin/time -v PISA parse2 -order -rule 'CB,R1:1-16;R1,R2' -1 reads_1.fq -q 0 atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L001_R2_001.fastq.gz,atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L002_R2_001.fastq.gz atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L001_R1_001.fastq.gz,atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L002_R1_001.fastq.gz -nw 2> parse_1.log 

/usr/bin/time -v PISA parse2 -order -rule 'CB,R1:1-16;R1,R2' -1 reads_2.fq -q 0 atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L001_R2_001.fastq.gz,atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L002_R2_001.fastq.gz atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L001_R3_001.fastq.gz,atac_pbmc_1k_nextgem_fastqs/atac_pbmc_1k_nextgem_S1_L002_R3_001.fastq.gz -nw 2> parse_2.log

# Alignment
bwa mem ../../pipeline/refdata-gex-GRCh38-2020-A/fasta/genome.fa reads_1.fq reads_2.fq -t 40 -o aln.sam

# Convert SAM to BAM
/usr/bin/time -v PISA sam2bam aln.sam -o aln.bam -@ 4 -report aln.csv 2> aln.csv

# Sort aligned BAM and index it
sambamba sort -t 4 aln.bam -o sorted.bam

# Convert to a fragment file
/usr/bin/time -v PISA bam2frag -o fragments.tsv.gz -cb CB -@ 4 sorted.bam 2>bam2frag.log

# Call peak with MACS2 in BED mode
macs2 callpeak -g hs --outdir ./ --format BED --call-summits --name pbmc --keep-dup all --nomodel --nolambda --shift -75 --extsize 150 --treatment fragments.tsv.gz
mkdir atac
/usr/bin/time -v PISA count2 -outdir atac -bed pbmc_peaks.narrowPeak -t 4 fragments.tsv.gz


########################
# Section Three
########################


# Test fastq convertion
cat 10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L*_R1*.gz > R1.fastq.gz
cat 10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L*_R2*.gz > R2.fastq.gz
/usr/bin/time -v umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --error-correct-cell --stdin  R1.fastq.gz --stdout R1_extract.fastq.g --read2-in R2.fastq.gz --read2-out R2_extract.fastq.gz --whitelist ./3M-february-2018.txt 1>umitools.log
/usr/bin/time -v PISA parse2 -rule 'CR,R1:1-16,3M-february-2018.txt.gz,CB,1;R1,R2' R1.fastq.gz R2.fastq.gz -1 reads_1.fq -t 1 -nw 1> parse_t1.log
/usr/bin/time -v PISA parse2 -rule 'CR,R1:1-16,3M-february-2018.txt,CB,1;R1,R2' R1.fastq.gz R2.fastq.gz -1 reads_1.fq -t 4 -nw 1> parse_t4.log

# Test gene annotation
wget -c https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam
/usr/bin/time -v ./Drop-seq_tools-2.5.1/TagReadWithGeneExonFunction I=10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam O=dropseq_anno.bam ANNOTATIONS_FILE=genes.gtf 2>dropseq.log
/usr/bin/time -v PISA anno -gtf genes.gtf 10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam -o anno1.bam -t 1 2>anno_t1.log
/usr/bin/time -v PISA anno -gtf genes.gtf 10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam -o anno1.bam -t 4 2>anno_t1.log
/usr/bin/time -v PISA anno -gtf genes.gtf 10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam -o anno1.bam -t 20 2>anno_t1.log
