#!/bin/bash
set -o pipefail
set -e
PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))
module load samtools
args=("$@")

function fragment_length {
	samtools view --threads $PROCESSORS $1 | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
	sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $2


};

fragment_length ${args[0]} ${args[1]}