#!/bin/bash

bamCompare -b1 samp.bam -b2 input.bam -o input_normalized.log2cpm.bw -p 6 --centerReads -bs 10 --blackListFileName blacklist.bed --normalizeUsing CPM --scaleFactorsMethod None