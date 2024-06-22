# WG_code_challenge
## Setup
For this challenge I used Python Anaconda for the environment and pip for dependency management. To create a new environment, execute the following in terminal:
```
conda create -n wg_code_challenge
```
Then activate the environment:
```
conda activate wg_code_challenge
```
and finally, install dependencies:
```
pip install -r requirements.txt
```

## Question 1
This is a coding task that involves finding the missing integer in a shuffled set of contiguous integers. The input contains a flat file that presents the series in the following format: 
There are multiple problems within a file. Individual problems are delimited with a header. The header is a greater than character $>$ followed by an integer $N$ as a header line. An example header is: $>10$. The value $N$ indicates the size of the set. Following this is a series of $N - 1$ integers which can span multiple lines.

### Solution
To run the solution for question one, execute the following command in terminal, from the project directory:
```
python question_one.py -f ABSOLUTE_FILE_PATH
```
where the `ABSOLUTE_FILE_PATH` is the absolute file path to the integer sequence input file. By default, the supplied file (`question_one.txt`) is called, with the following:
```
python question_one.py
```

Output results will have the following format:
```
Missing element: ELEMENT_NUMBER
...
```
for each sequence in the file.

### Testing
To run the test suite, execute the following command in terminal, from the project directory:
```
python test_question_one.py
```

## Question 2
### Requirement
Find the region in the toy genome with 2x copy number gain.

**Assets:**
* toy_genome.fa: The reference genome.
* toy_genome_genes.bed: Annotations indicating mock genes.
* baseline_reads_R*.fastq: Simulated reads, unmodified.
* reads_R*.fastq: Simulated reads, modified to contain a 2x copy gain somewhere.

### Assumptions
1. All fastq files contain a quality rating of `2` for all nucleotide readings. Therefore, quality can be ignored in this analysis.
2. Readings are not amplified and contain no noise:
    * Readings have not been amplied with Polymerase Chain Reactions.
    * Readings do not contain optical duplications.

### Strategy 1
The steps to accomplish the task are as follows:
1. Use existing Burrows-Wheeler implementation (bwa) to align reads with reference genome.
2. Use Sequence Alignment Map tools (samtools) to determine depth of coverage.
3. Normalize coverage distributions so that samples can be reliably compared.
4. Determine copy number for coverage distributions.
5. Apply bounded filter to remove regions not of interest.
6. Examine positions of gains within region interest.
7. Correlate largest contiguous region with gene annotations to find most relevant gene. This will most likely be the gene where the copy gain occurred.

Steps 1 and 2 were executed using command line tools in WSL Ubuntu.

#### Step 1 - Align Reads
First, I installed the bwa command line tool.
```shell
# Update references
sudo apt update

# Install BWA
sudo apt install bwa
```

From this I indexed the toy genome so that reads can be appropriately aligned. Then the baseline and modified reads were aligned separately.

```shell
# Create Toy Genome Indices
bwa index toy_genome.fa

# Create Baseline Reads Alignment
bwa mem toy_genome.fa baseline_reads_R1.fastq baseline_reads_R2.fastq | gzip -3 > baseline_alignment.sam.gz

# Create Modified Reads Alignment
bwa mem toy_genome.fa reads_R1.fastq reads_R2.fastq | gzip -3 > modified_alignment.sam.gz
```

#### Step 2 - Determine Coverage Depth
First, I installed samtools.

```shell
# Install Samtools
sudo apt install samtools
```

Then I converted the alignments to a binary format, which was necessary for sorting.

```shell
# Convert SAM to BAM
samtools view -S -b baseline_alignment.sam.gz > baseline_alignment.bam.gz
samtools view -S -b modified_alignment.sam.gz > modified_alignment.bam.gz

# Sort BAM files
samtools sort baseline_alignment.bam.gz -o sorted_baseline_alignment.bam.gz
samtools sort modified_alignment.bam.gz -o sorted_modified_alignment.bam.gz
```

Finally, I calculated the coverage depth for each alignment separately.

```shell
# Determine depth of coverage
samtools depth -a sorted_baseline_alignment.bam.gz > baseline_coverage.txt
samtools depth -a sorted_modified_alignment.bam.gz > modified_coverage.txt
```

Steps 3-7 of the solution continue in the `Question 2 - Strategy 1` jupyter notebook.

**NOTE:** All files generated in the solution pipeline are included in the `data` directory as proof that the solution was executed as documented.

### Strategy 2
Create custom tools to align reads and compare maps to determine where 2x copy gain is located.