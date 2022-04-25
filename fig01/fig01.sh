# Download test dataset
wget -c https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam
PISA rmdup -tags CB,UB -q 20 -nw Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam  -o rmdup.bam -@ 10
sambamba index rmdup.bam

bcftools mpileup -Ou -f ../../pipeline/refdata-gex-GRCh38-2020-A/fasta/genome.fa rmdup.bam | bcftools call -vmO z -o var.vcf.gz

# Generate exon regions from GTF
awk '{if($3 == "exon"){printf("%s\t%d\t%d\t.\t.\t%c\n",$1,$4,$5,$7)}}' ../../pipeline/refdata-gex-GRCh38-2020-A/genes/genes.gtf | sort -k1,1 -k2,2n -k3,3n |uniq |bedtools merge -s -c 6 -o collapse > tmp.bed
awk '{printf("%s\t%d\t%d\t.\t.\t%c\n",$1,$2,$3,$4)}' tmp.bed > isoforms.bed

# Download RepeatMask result from UCSC
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
# Convert to BED
zcat rmsk.txt.gz|awk -vFS="\t" '{printf("%s\t%d\t%d\t%s\t.\t%c\t%s\n",$6,$7,$8,$11,$10,$12)}' > rmsk.bed
# Select TEs
awk 'BEGIN{n=split("DNA LINE LTR SINE Satellite Retroposon", TE); for (i in TE) {val[TE[i]]=""}}; ($7 in val) {print $0}' rmsk.bed > TE.bed

mkdir exp
mkdir var
mkdir isoform
mkdir TE

# Feature annotation
PISA anno -bed isoforms.bed -tag IS rmdup.bam -o anno_isoform.bam -t 10 
PISA anno -bed TE.bed rmdup.bam -o anno_te.bam -tag TE -t 10
PISA anno -vcf ./var.vcf.gz -vtag VF rmdup.bam -o anno_vars.bam -t 10

# Count matrix
PISA count -cb CB -anno-tag GN -outdir exp rmdup -@ 10 -umi UB 
PISA count -cb CB -anno-tag RM -outdir TE anno_te.bam -@ 10 -umi UB 
PISA count -cb CB -anno-tag IS -outdir isoform anno_isoform.bam -@ 10 -umi UB 
PISA count -cb CB -anno-tag VF -outdir var anno_vars.bam -@ 10 -umi UB 
