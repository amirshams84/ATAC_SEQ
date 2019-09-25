#! /bin/bash
set -o pipefail
set -e

module load samtools
module load bedtools
args=("$@")
function frip_score {
	printf "%s\t%s\t%s\t%s\n" "Sample_name" "total_read" "total_peak" "read_in_peak_count" > $4
	read_in_peak_count=$(bedtools intersect -a <(zcat -f $2) -b <(zcat -f $3) -wa -u | wc -l)
	total_read=$(zcat $2 | wc -l )
	total_peak=$(zcat $3 | wc -l )
	printf "%s\t%s\t%s\t%s\n" "$1" "$total_read" "$total_peak" "$read_in_peak_count" >> $4

};


frip_score ${args[0]} ${args[1]} ${args[2]} ${args[3]}