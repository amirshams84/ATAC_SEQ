#! /bin/bash
set -o pipefail
set -e
args=("$@")
if [ $# -eq 0 ]; then
    echo "No JSON file provided!!!"
    echo "Aborting!!"
    exit 1
fi


function read_link() {
    local path=$1
    if [ -d $path ] ; then
        local abspath=$(cd $path; pwd -P)
    else
        local prefix=$(cd $(dirname -- $path) ; pwd -P)
        local suffix=$(basename $path)
        local abspath="$prefix/$suffix"
    fi
    if [ -e $abspath ] ; then
        echo $abspath
    else
        echo "$1 is not accessible"
		echo "make sure the path is correct(Please use absoulte path)!"
		echo "Aborting!!"
    	exit 1
    fi
};

echo "Pipeline execution initiated at: "$(date)

JSON_CONFIG_FILE=${args[0]}
JSON_CONFIG_FILE=$(read_link $JSON_CONFIG_FILE)
####################################
#
EXECUTION_MODE="DEVELOPMENT"
EXECUTION_MODE="DEPLOYMENT"
EXECUTION_MODE="TEST"


EXECUTION_SCRIPT=/data/shamsaddinisha/ATAC_Seq/Snakemake/execution/atac_seq_execution.sh

bash /data/shamsaddinisha/ATAC_Seq/Snakemake/execution/atac_seq_execution.sh $JSON_CONFIG_FILE $EXECUTION_MODE

echo "Pipeline execution successfully finished at: "$(date)



