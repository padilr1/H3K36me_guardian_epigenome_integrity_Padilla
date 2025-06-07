#!/bin/bash

# after aligning and processing the FASTQ file, depth-normalized bigWigs (in counts per million (CPM)) for each sample was generated using:
bamCoverage -b ${sample}.sorted.bam -o ${sample}.cpm.bw --normalizeUsing CPM --centerReads -e 200 -bs 10

# input-normalized bigWigs were generated using:
bamCompare -b1 ${ChIP_sample}.bam -b2 ${input_sample}.bam -o ${sample}.input_normalized.log2cpm.bw -e 200 --centerReads -bs 10 --blackListFileName blacklist.bed --normalizeUsing CPM --scaleFactorsMethod None

# samples were merged using:
bigwigCompare --operation mean -bs 10 -b1 ${sample_rep1}.cpm.bw -b2 ${sample_rep2}.cpm.bw -o merged.bw
