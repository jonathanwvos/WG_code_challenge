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
Steps:
1. Use existing Burrows-Wheeler implementation to align reads with reference genome.
2. Use Sequence Alignment Map tools to determine depth of coverage.
3. Normalize coverage distributions so that samples can be reliably compared.
4. Determine copy number for coverage distributions.
5. Apply bounded filter to remove regions not of interest.
6. Examine positions of gains within region interest.
7. Correlate largest contiguous region with gene annotations to find most relevant gene. This will most likely be the gene where the copy gain occurred.

Refer to the `Question 2 - Strategy 1` notebook for the implementation.

### Strategy 2
Create custom tools to align reads and compare maps to determine where 2x copy gain is located.