#! /bin/bash
set -o pipefail
set -e
args=("$@")
if [ $# -eq 0 ]; then
    echo "No Input file provided!!!"
    echo "Aborting!!"
    exit 1

fi
####################################################
#

declare -a json_List
function parse_json {
	
	index=0
	IFS="&&"
	export PYTHONIOENCODING=utf8
	for line in $(cat $1 | python -c 'import os,sys,json; \
		data = json.load(sys.stdin); \
		print(\
		os.path.expanduser(data["DIRECTORY"]["WORKDIR"])+"&&"+\
		str(data["GENERAL"]["PROJECT"])+"&&"+\
		str(data["GENERAL"]["EXPERIMENT"])+"&&"+\
		str(data["GENERAL"]["TITLE"])
		)')
	do
		#echo $index
		#echo $line
		json_List[$index]=$line
		index=$(($index+1))
		
	done
	IFS=" "
	
};

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
####################################################
#
#Parse json file
JSON_CONFIG_FILE=${args[0]}
#
parse_json $JSON_CONFIG_FILE
WORK_DIR=${json_List[0]}
mkdir -p $WORK_DIR
WORK_DIR=$(read_link $WORK_DIR)
#
PROJECT=${json_List[2]}
EXPERIMENT=${json_List[4]}
TITLE=${json_List[6]}
#
LOG_DIR=$WORK_DIR/${PROJECT}/${EXPERIMENT}/${TITLE}/logs
mkdir -p $LOG_DIR
LOG_DIR=$(read_link $LOG_DIR)
#
ATAC_SEQ_DIR=/data/shamsaddinisha/ATAC_Seq
SNAKEMAKE_DIR=${ATAC_SEQ_DIR}/Snakemake
SNAKE_FILE_DIR=$SNAKEMAKE_DIR/snakefile
CLUSTER_CONFIG=${ATAC_SEQ_DIR}/Snakemake/cluster_config
####################################################
#

EXECUTION_MODE=${args[1]}

if [ "$EXECUTION_MODE" == "DEVELOPMENT" ]
then
	PROCESSORS=10
	MEMORY=10000
elif [ "$EXECUTION_MODE" == "TEST" ]
then
	PROCESSORS=10
	MEMORY=50000
	CLUSTER_CONFIG_FILE=${CLUSTER_CONFIG}/cluster_test.yaml
elif [ "$EXECUTION_MODE" == "DEPLOYMENT" ]
then
	PROCESSORS=30
	MEMORY=100000
	CLUSTER_CONFIG_FILE=${CLUSTER_CONFIG}/cluster_deploy.yaml
fi
####################################################
#

module load snakemake || exit 1
####################################################
#
snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --unlock

if [ "$EXECUTION_MODE" == "DEVELOPMENT" ]
then

	snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/alignment.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/post_alignment.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/peak_calling.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/peak_overlap.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/peak_annotate.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

	snakemake --snakefile $SNAKE_FILE_DIR/nucleoatac.py --configfile $JSON_CONFIG_FILE --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete

else

	mkdir -p ${LOG_DIR}/pre_process
	snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/pre_process/{cluster.output} --error=${LOG_DIR}/pre_process/{cluster.error} {cluster.extra}"
	
	mkdir -p ${LOG_DIR}/alignment
	snakemake --snakefile $SNAKE_FILE_DIR/alignment.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/alignment/{cluster.output} --error=${LOG_DIR}/alignment/{cluster.error} {cluster.extra}"

	mkdir -p ${LOG_DIR}/post_alignment
	snakemake --snakefile $SNAKE_FILE_DIR/post_alignment.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/post_alignment/{cluster.output} --error=${LOG_DIR}/post_alignment/{cluster.error} {cluster.extra}"

	mkdir -p ${LOG_DIR}/peak_calling
	snakemake --snakefile $SNAKE_FILE_DIR/peak_calling.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/peak_calling/{cluster.output} --error=${LOG_DIR}/peak_calling/{cluster.error} {cluster.extra}"

	mkdir -p ${LOG_DIR}/peak_overlap
	snakemake --snakefile $SNAKE_FILE_DIR/peak_overlap.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/peak_overlap/{cluster.output} --error=${LOG_DIR}/peak_overlap/{cluster.error} {cluster.extra}"

	mkdir -p ${LOG_DIR}/peak_annotate
	snakemake --snakefile $SNAKE_FILE_DIR/peak_annotate.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/peak_annotate/{cluster.output} --error=${LOG_DIR}/peak_annotate/{cluster.error} {cluster.extra}"

	mkdir -p ${LOG_DIR}/nucleoatac
	snakemake --snakefile $SNAKE_FILE_DIR/nucleoatac.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/nucleoatac/{cluster.output} --error=${LOG_DIR}/nucleoatac/{cluster.error} {cluster.extra}"

fi
#
####################################################