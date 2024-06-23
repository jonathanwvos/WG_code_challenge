# Alignment Statistics
Using the `samtools` command line utility, we can draft statistics of the baseline and modified alignment steps. Doing so will endow us with confidence that subsequent steps in the analyis are reliable.

## Baseline Alignment
The following command was used to generate read alignemnt stats for baseline reads.

```shell
samtools flagstat baseline_alignment.bam.gz
```

Which yielded the following results:

```
199932 + 0 in total (QC-passed reads + QC-failed reads)
199932 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
199932 + 0 mapped (100.00% : N/A)
199932 + 0 primary mapped (100.00% : N/A)
199932 + 0 paired in sequencing
99966 + 0 read1
99966 + 0 read2
199692 + 0 properly paired (99.88% : N/A)
199932 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

It seems clear that the alignment was impecable, with no reads failing QC checks.

## Modified Alignment

Similarly, for the the modified reads:

```shell
samtools flagstat modified_alignment.bam.gz
```

Yielding similar, unsurprising results!

```
201094 + 0 in total (QC-passed reads + QC-failed reads)
201094 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
201092 + 0 mapped (100.00% : N/A)
201092 + 0 primary mapped (100.00% : N/A)
201094 + 0 paired in sequencing
100547 + 0 read1
100547 + 0 read2
200840 + 0 properly paired (99.87% : N/A)
201092 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

It is safe to say that the alignment can be properly relied on in subsequent steps in the analysis.

# Mapping Quality
In addition to alignment quality, we can also inspect the mapping quality of reads to reference genome. Once again, we use `samtools` to assist in the analysis.

The following reference guide can be used to interpret MAPQ scores:
* High MAPQ Scores (e.g., 30-60): Indicates high confidence in the alignment.
* Moderate MAPQ Scores (e.g., 10-29): Indicates moderate confidence. Reads with these scores might require further inspection or filtering.
* Low MAPQ Scores (e.g., 0-9): Indicates low confidence. These reads are often ambiguously aligned or have multiple equally likely alignment positions.

## Baseline Alignment
Generate a distribution of MAPQ scores for the baseline reads alignment:

```shell
samtools view baseline_reads.bam.gz | awk '{print $5}' | sort | uniq -c > baseline_mapq_distribution.txt
```

The following was copied from the `baseline_mapq_distribution.txt` file:

```
    2 54
    3 58
199927 60
```

All reads have high scores, with only 5 scoring below 60.

## Modified Alignment
Similarly, for the modified reads alignment:

```shell
samtools view modified_alignment.bam.gz | awk '{print $5}' | sort | uniq -c > modified_mapq_distribution.txt
```

The following was copied from the `modified_mapq_distribution.txt` file:

```
    2 0
    1 42
    3 58
201088 60
```
The vast majority of reads have high scores, with only 4 scoring just below 60, and 2 scoring 0.
