#!/bin/bash
#SBATCH --job-name=STARmap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/HeLa_out.%j
#SBATCH --error=/data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/HeLa_err.%j
echo "HeLa processing start time: $(date)";
echo ""
STAR --runThreadN 16 --genomeDir  /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/STAR/hg38 --genomeLoad LoadAndExit --outSAMtype None
echo ".......... loaded genome: $(date)";
echo ""
echo "Mapping G1 sample of 2018_02_09_ILL_DGB project to hg38"
echo ""
echo ".......... G1 mapping start time: $(date)"
STAR --readFilesIn /data/nieduszynski/NGS/2018_02_09_ILL_DGB/6-G1-HeLa_S2_L001_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/6-G1-HeLa_S2_L002_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/6-G1-HeLa_S2_L003_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/6-G1-HeLa_S2_L004_R1_001.fastq.gz \
--runThreadN 16 \
--alignIntronMax 1 \
--genomeDir  /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/STAR/hg38 \
--genomeLoad LoadAndKeep \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1_ \
--readFilesCommand gunzip -c -k \
--alignEndsType EndToEnd \
--outFilterScoreMin 30 \
--outFilterMismatchNmax 3 \
--outFilterMultimapScoreRange 251 \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate
echo ".......... G1 mapping end time: $(date)";
echo ""
echo "Mapping S sample of 2018_02_09_ILL_DGB project to hg38"
echo ""
echo ".......... S mapping start time: $(date)"
echo ""
STAR --readFilesIn /data/nieduszynski/NGS/2018_02_09_ILL_DGB/12-S-HeLa_S1_L001_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/12-S-HeLa_S1_L002_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/12-S-HeLa_S1_L003_R1_001.fastq.gz,\
/data/nieduszynski/NGS/2018_02_09_ILL_DGB/12-S-HeLa_S1_L004_R1_001.fastq.gz \
--runThreadN 16 \
--alignIntronMax 1 \
--genomeDir  /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/STAR/hg38 \
--genomeLoad LoadAndKeep \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/G1_ \
--readFilesCommand gunzip -c -k \
--alignEndsType EndToEnd \
--outFilterScoreMin 30 \
--outFilterMismatchNmax 3 \
--outFilterMultimapScoreRange 251 \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate
echo ".......... S mapping end time: $(date)";
echo ""
STAR --runThreadN 16 --genomeDir  /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/STAR/hg38 --genomeLoad Remove
echo ".......... unloaded genome: $(date)";
echo ""
samtools index /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1_Aligned.sortedByCoord.out.bam
samtools index /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S_Aligned.sortedByCoord.out.bam
test ! -e /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed && samtools idxstats /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1_Aligned.sortedByCoord.out.bam | awk 'BEGIN {OFS="\t"} {if ($2>0) print ($1,$2)}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/hg38.txt && bedtools makewindows -g /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/hg38.txt -w 10000 > /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed
(samtools view -H /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1_Aligned.sortedByCoord.out.bam; samtools view -h -@ 15 -q 30 -F 3840 -L /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1_Aligned.sortedByCoord.out.bam | grep -w "NH:i:1") | samtools view -@ 15 -b - | bedtools genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS="\t"} {if ($3>0) print $1,$2,$2,"name",$3}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1-HeLa_coverage.bed
bedtools map -a /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed -b /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1-HeLa_coverage.bed -null 0 -o sum | awk 'BEGIN {OFS="\t"} {if ($4>0) print $1,$2,$3,"name",$4}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/processed/G1-HeLa_10kb.bed
(samtools view -H /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S_Aligned.sortedByCoord.out.bam; samtools view -h -@ 15 -q 30 -F 3840 -L /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S_Aligned.sortedByCoord.out.bam | grep -w "NH:i:1") | samtools view -@ 15 -b - | bedtools genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS="\t"} {if ($3>0) print $1,$2,$2,"name",$3}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S-HeLa_coverage.bed
bedtools map -a /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.10000bp.bed -b /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S-HeLa_coverage.bed -null 0 -o sum | awk 'BEGIN {OFS="\t"} {if ($4>0) print $1,$2,$3,"name",$4}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/processed/S-HeLa_10kb.bed
bedtools map -a /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.50000bp.bed -b /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/tmp/G1-HeLa_coverage.bed -null 0 -o sum | awk 'BEGIN {OFS="\t"} {if ($4>0) print $1,$2,$3,"name",$4}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/G1-HeLa/processed/G1-HeLa_50kb.bed
bedtools map -a /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/genomeWindows/hg38.50000bp.bed -b /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/tmp/S-HeLa_coverage.bed -null 0 -o sum | awk 'BEGIN {OFS="\t"} {if ($4>0) print $1,$2,$3,"name",$4}' > /data/nieduszynski/NGS/2018_02_09_ILL_DGB/STAR/S-HeLa/processed/S-HeLa_50kb.bed
echo ".......... HeLa processing end time: $(date)";
echo ""
