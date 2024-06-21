# Update references
sudo apt update

# Install BWA
sudo apt install bwa

# Create Toy Genome Indices
bwa index toy_genome.fa

# Create Baseline Reads Alignment
bwa mem toy_genome.fa baseline_reads_R1.fastq baseline_reads_R2.fastq | gzip -3 > baseline_alignment.sam.gz

# Create Modified Reads Alignment
bwa mem toy_genome.fa reads_R1.fastq reads_R2.fastq | gzip -3 > modified_alignment.sam.gz

# Install Samtools
sudo apt install samtools

# Convert SAM to BAM
samtools view -S -b baseline_alignment.sam.gz > baseline_alignment.bam.gz
samtools view -S -b modified_alignment.sam.gz > modified_alignment.bam.gz

# Sort BAM files
samtools sort baseline_alignment.bam.gz -o sorted_baseline_alignment.bam.gz
samtools sort modified_alignment.bam.gz -o sorted_modified_alignment.bam.gz

# Determine depth of coverage
samtools depth -a sorted_baseline_alignment.bam.gz > baseline_coverage.txt
samtools depth -a sorted_modified_alignment.bam.gz > modified_coverage.txt