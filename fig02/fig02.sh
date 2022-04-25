# Download test data
samtools view https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam chr11:18227869-18270672 -b -o SAA1-4.bam

# Convert to a FASTQ+ file
PISA bam2fq -tags CB,UB -f SAA1-4.bam -o pick.fq

# Order FASTQ+ in cell blocks
PISA fsort -tags CB pick.fq -o sort.fq.gz

# Assembly each block with Trinity
PISA stream -tags CB -script 'Trinity --seqType fq --SS_lib_type F --single ${FQ} --max_memory 1G 2>/dev/null 1>/dev/null; seqtk rename trinity_out_dir.Trinity.fasta ${UBI}_ 2>/dev/null' -t 10 -fa  -nw ./sort.fq.gz -o asm.fa

# Remap assembled reads to reference
minimap2 -x splice -a ~/Documents/datasets/GRCh38/fasta/genome.fasta asm.fa 1> asm_aln.sam

# Convert SAM to BAM
PISA sam2bam asm_aln.sam -o asm_aln.bam

# Sort for IGV
sambamba sort asm_aln.bam -o asm_sorted.bam

