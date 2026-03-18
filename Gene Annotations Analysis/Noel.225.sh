# Dr. Singh Section
#1a
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gtf.gz > B73_maize.gtf.gz
gunzip B73_maize.gtf.gz > B73_maize.gtf
# final output is the unzipped B73_maize.gtf.gz file that was downloaded from the link.

#1b
awk -F'\t' '$3=="gene"{if(match($9,/gene_biotype "([^"]+)"/,a)) print a[1]; else print "NA"}' B73_maize.gtf | sort | uniq -c | sort -nr
# final output for maize GTF is:
#  34327 protein_coding
#   5348 lncRNA
#   3887 pseudogene
#   2439 rRNA
#   1366 tRNA
#   1323 transcribed_pseudogene
#    655 snoRNA
#    241 snRNA
#    164 miRNA
#    139 misc_RNA
#      8 tRNA_pseudogene
#      4 ncRNA_pseudogene
#      1 antisense_RNA

curl -s https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz > Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz > Homo_sapiens.GRCh38.115.gtf
awk -F'\t' '$3=="gene"{if(match($9,/gene_biotype "([^"]+)"/,a)) print a[1]; else print "NA"}' Homo_sapiens.GRCh38.115.gtf | sort | uniq -c | sort -nr
#  final output for human GTF is:
#  35034 lncRNA
#  20121 protein_coding
#   9489 processed_pseudogene
#   2216 misc_RNA
#   1949 unprocessed_pseudogene
#   1910 snRNA
#   1879 miRNA
#   1587 transcribed_unprocessed_pseudogene
#   1149 transcribed_processed_pseudogene
#   1019 TEC
#    942 snoRNA
#    497 rRNA_pseudogene
#    201 transcribed_unitary_pseudogene
#    187 IG_V_pseudogene
#    146 IG_V_gene
#    107 TR_V_gene
#     89 unitary_pseudogene
#     79 TR_J_gene
#     53 rRNA
#     49 scaRNA
#     37 IG_D_gene
#     33 TR_V_pseudogene
#     22 Mt_tRNA
#     19 artifact
#     18 IG_J_gene
#     14 IG_C_gene
#      9 IG_C_pseudogene
#      8 ribozyme
#      6 TR_C_gene
#      5 TR_D_gene
#      5 sRNA
#      4 vault_RNA
#      4 TR_J_pseudogene
#      4 pseudogene
#      3 IG_J_pseudogene
#      2 translated_processed_pseudogene
#      2 Mt_rRNA
#      1 IG_pseudogene

#1e
awk -F'\t' '$1=="NC_050105.1" && $3=="gene" && match($9,/gene_biotype "([^"]+)"/,b) && b[1]=="protein_coding" && match($9,/gene_id "([^"]+)"/,g){print $1"\t"($4-1)"\t"$5"\t"g[1]"\t.\t"$7}' B73_maize.gtf > NC_050105.1_protein_coding_genes.bed
# final output is the 6-column BED file lists the chromosome NC_050105.1 in the first column, the start position of that gene's coordinates in the second column, the end position of that gene's coordinates in the third column, the gene_id in the fourth column, a placeholder score (either 0 or .) is in the fifth column, the strand notation (+ or -) is in the sixth column.

#2a
sinteractive -A PAS3124 -c 48 -t 1:00:00
module load fastqc/0.12.1
mkdir -p fastqc_out
# final output is a directory for the fastqc analyses ran.

fastqc -t 8 -o fastqc_out c1_1.fastq c1_2.fastq c2_1.fastq c2_2.fastq m1_1.fastq m1_2.fastq m2_1.fastq m2_2.fastq
# final outputs are the .html files and compressed .zip files for the analyses.

#2b
module load trimmomatic/0.38
trimmomatic PE c1_1.fastq c1_2.fastq c1_paired_F.fastq c1_unpaired_F.fastq c1_paired_R.fastq c1_unpaired_R.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
# final output is the forward and reverse, paired and unpaired, trimmed files for the first control (c1).

trimmomatic PE c2_1.fastq c2_2.fastq c2_paired_F.fastq c2_unpaired_F.fastq c2_paired_R.fastq c2_unpaired_R.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
# final output is the forward and reverse, paired and unpaired, trimmed files for the second control (c2).

trimmomatic PE m1_1.fastq m1_2.fastq m1_paired_F.fastq m1_unpaired_F.fastq m1_paired_R.fastq m1_unpaired_R.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
# final output is the forward and reverse, paired and unpaired, trimmed files for the first mutant (m2).

trimmomatic PE m2_1.fastq m2_2.fastq m2_paired_F.fastq m2_unpaired_F.fastq m2_paired_R.fastq m2_unpaired_R.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
# final output is the forward and reverse, paired and unpaired, trimmed files for the second mutant (m2).

# Dr. Ou Section
#1a
sinteractive -A PAS3124 -n 48 -t 1:00:00
# Starts interactive node of 48 CPUs at OSC

# set up the base environment
module load miniconda3/24.1.2-py310
conda env create -f helixer.yml
conda activate helixer

# install helixer
python -m pip install helixerlite

# test if you have helixerlite installed
helixerlite -h

#Annotate Maize sequence with land_plant model
module load miniconda3/24.1.2-py310
conda activate helixer
nohup helixerlite --fasta B73_chr3_10M.fasta --lineage land_plant --out B73_chr3_10M.land_plant.output.gff3 -c 30 &  

#1b
grep -P '\s+gene\s+' B73_chr3_10M.land_plant.output.gff3 > helixer_genes.gff3
grep -P '\s+gene\s+' B73_NC_050098.1_10M.Zm00001eb.1.genes.gff3 > ref_genes.gff3

#1c 
# GFF conversion to BED
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9, ".", $7}' helixer_genes.gff3 > helixer_genes.bed
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9, ".", $7}' ref_genes.gff3     > ref_genes.bed

# install BEDtools
conda install -c bioconda bedtools

# Avoiding naming mismatches
awk '$1=="NC_050098.1"' helixer_genes.bed > helixer_chr3.bed
awk '$1=="NC_050098.1"' ref_genes.bed     > ref_chr3.bed

# Finding overlaps are between the reference and Helixer annotation
bedtools intersect -a helixer_chr3.bed -b ref_chr3.bed -f 0.90 -F 0.90 -wa -wb | wc -l

#2a
# install EDTA
git clone https://github.com/oushujun/EDTA.git  
module load miniconda3/24.1.2-py310

conda create -n EDTA -y  
conda activate EDTA  
conda install -c conda-forge -c bioconda edta -y

conda install -c bioconda rmblast repeatmasker -y

perl ./EDTA/EDTA.pl --genome ./B73_chr3_10M.fasta --curatedlib ./EDTA/database/maizeTE11122019 --anno 1 --threads 48

#2b
B73_chr3_10M.fasta.mod.EDTA.TEanno.sum

#2c
cat B73_chr3_10M.fasta.mod.EDTA.TEanno.gff3 | grep 'method=structural' | less

#3a
sinteractive -A PAS3124 -c 48 -t 1:00:00
cd ~/MG5645/noel.225/MG5795-2025
git pull
module load miniconda3/24.1.2-py310
conda create -n sv -c bioconda -c conda-forge survivor "sniffles>2" assemblytics minimap2 bcftools samtools bedtools   
conda activate sv

nucmer -maxmatch -l 100 -c 500 B73_chr3_10M.fasta Mo17_chr3_10M.wrap50.fasta

Assemblytics out.delta B73-Mo17_chr3_10M 100 50 50000

less B73-Mo17_chr3_10M.Assemblytics_structural_variants.bed
less B73-Mo17_chr3_10M.Assemblytics_structural_variants.summary

SURVIVOR convertAssemblytics B73-Mo17_chr3_10M.Assemblytics_structural_variants.bed 10 B73_chr3_10M.WGA.SV.vcf
less -S B73_chr3_10M.WGA.SV.vcf

#3b
grep -v '^#' B73_chr3_10M.WGA.SV.vcf | wc -l

#3c
grep -v '^#' B73_chr3_10M.WGA.SV.vcf | \
awk 'BEGIN{OFS="\t"}{
  chrom=$1
  start=$2-1
  if (match($8,/END=([0-9]+)/,a)) end=a[1]
  else end=$2
  if (match($8,/SVTYPE=([^;]+)/,b)) svtype=b[1]
  else svtype="SV"
  print chrom,start,end,svtype
}' > sv_chr3.bed

awk '$1=="NC_050098.1"' \
  B73_chr3_10M.fasta.mod.EDTA.TEanno.gff3 | \
awk 'BEGIN{OFS="\t"}{
  print $1,$4-1,$5,$9
}' > te_chr3.bed

bedtools intersect \
  -a sv_chr3.bed \
  -b te_chr3.bed \
  -f 0.90 -F 0.90 -wa -wb | wc -l
  
  bedtools intersect \
  -a sv_chr3.bed \
  -b te_chr3.bed \
  -f 0.90 -F 0.90 -wa -wb | \
cut -f4 | sort | uniq -c

#4a
# Map Mo17 long reads to the B73 reference sequence
minimap2 -ax map-hifi -t 48 -R '@RG\tID:Mo17\tSM:Mo17\tLB:Mo17\tPL:PACBIO' B73_chr3_10M.fasta Mo17.SRR30894834.chr3_10M.hifi_reads.fastq.gz | samtools sort -@ 8 -o B73_chr3_10M.Mo17.sorted.bam

# index  bam file 
samtools index B73_chr3_10M.Mo17.sorted.bam

# SV calling
sniffles \
  --input B73_chr3_10M.Mo17.sorted.bam \
  --vcf B73_chr3_10M.Mo17.SV.vcf \
  --reference B73_chr3_10M.fasta \
  --minsvlen 10 \
  --minsupport 3 \
  --mapq 40
	
# Inspecting Files
less -S B73_chr3_10M.Mo17.SV.vcf
grep -vc '#' B73_chr3_10M.Mo17.SV.vcf

#4b
bcftools filter -i 'QUAL>=50 && ABS(INFO/SVLEN)<=20000 && GT="1/1"' B73_chr3_10M.Mo17.SV.vcf > B73_chr3_10M.Mo17.filtered.SV.vcf
grep -vc '#' B73_chr3_10M.Mo17.filtered.SV.vcf

#5a
sinteractive -A PAS3124 -c 48 -t 1:00:00
cd ~/MG5645/noel.225/MG5795-2025
git pull

module load miniconda3/24.1.2-py310
conda activate mapping
conda install pbmm2 pb-cpg-tools -c bioconda -y

samtools view B73.SRR11606869.chr3_10M.hifi_reads.bam | less

pbmm2 index B73_chr3_10M.fasta B73_chr3_10M.mmi  
pbmm2 align B73_chr3_10M.fasta B73.SRR11606869.chr3_10M.hifi_reads.bam B73_chr3_10M.hifi_reads.aln.bam --sort --num-threads 30 -J 2 --preset CCS

aligned_bam_to_cpg_scores --bam B73_chr3_10M.hifi_reads.aln.bam \
	--modsites-mode reference --ref B73_chr3_10M.fasta \
	--pileup-mode model --min-mapq 60 --min-coverage 4 \
	--output-prefix B73_chr3_10M.hifi_reads.CpG --threads 30
	
#5b
awk 'BEGIN{OFS="\t"} NF>=9 && $0 !~ /^#/ && tolower($0) !~ /helitron/ {print $1,$4-1,$5,$9}' B73_chr3_10M.fasta.mod.EDTA.TEanno.gff3 > B73_chr3_10M.TE.noHelitron.bed
  
bedtools intersect  -a B73_chr3_10M.hifi_reads.CpG.combined.UMR.bed -b B73_chr3_10M.TE.noHelitron.bed -f 1.0 -wa -wb | awk '{print $NF}' | sort -u | wc -l
# final output of 23

#5c
awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="gene" {print $1,$4-1,$5,$9,".",$7}' B73_chr3_10M.land_plant.output.gff3 > helixer_genes.bed
  
bedtools intersect -a B73_chr3_10M.hifi_reads.CpG.combined.UMR.bed -b helixer_genes.bed -f 1.0 -wa -wb | awk '{print $NF}' | sort -u | wc -l
