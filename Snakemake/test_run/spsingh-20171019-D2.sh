#! /bin/bash
set -o pipefail
set -e
echo "Pipeline execution initiated at: "$(date)
##############################################

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))
##############################################

WORK_DIR=/data/shamsaddinisha/jfarber-20160119-P1/jfarber-20160119-E1
JSON_CONFIG_FILE=/data/shamsaddinisha/jfarber-20160119-P1/jfarber-20160119-E1/script/spsingh-20171019-D2.json
LOG_DIR=$WORK_DIR/logs
SNAKEMAKE_DIR=/data/shamsaddinisha/ATAC_Seq/Snakemake
CLUSTER_CONFIG_DIR=$SNAKEMAKE_DIR/cluster_config
SNAKE_FILE_DIR=$SNAKEMAKE_DIR/snakefile
####################################################

mkdir -p $LOG_DIR/pre_process
mkdir -p $LOG_DIR/alignment
mkdir -p $LOG_DIR/post_alignment
mkdir -p $LOG_DIR/peak_calling
mkdir -p $LOG_DIR/peak_overlap
mkdir -p $LOG_DIR/peak_annotate
mkdir -p $LOG_DIR/peak_analysis
#mkdir -p ./logs/build_trackhub
########################################

module load snakemake || exit 1
module load java

snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --unlock
#
snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/pre_process.yaml \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/pre_process/{cluster.output} --error=$LOG_DIR/pre_process/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#
snakemake --snakefile $SNAKE_FILE_DIR/alignment.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/alignment.yaml  \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/alignment/{cluster.output} --error=$LOG_DIR/alignment/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#

snakemake --snakefile $SNAKE_FILE_DIR/post_alignment.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/post_alignment.yaml \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/post_alignment/{cluster.output} --error=$LOG_DIR/post_alignment/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#
snakemake --snakefile $SNAKE_FILE_DIR/peak_calling.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/peak_calling.yaml \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/peak_calling/{cluster.output} --error=$LOG_DIR/peak_calling/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#
snakemake --snakefile $SNAKE_FILE_DIR/peak_overlap.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/peak_overlap.yaml \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/peak_overlap/{cluster.output} --error=$LOG_DIR/peak_overlap/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#
snakemake --snakefile $SNAKE_FILE_DIR/peak_annotate.py --configfile $JSON_CONFIG_FILE  --local-cores=$PROCESSORS --cores=$PROCESSORS --cluster-config $CLUSTER_CONFIG_DIR/peak_annotate.yaml \
-j 100 --latency-wait 120 --keep-going --rerun-incomplete --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
--job-name={cluster.jobname} --output=$LOG_DIR/peak_annotate/{cluster.output} --error=$LOG_DIR/peak_annotate/{cluster.error} --time={cluster.time} --mail-type=FAIL {cluster.extra}"
#
#snakemake --snakefile peak_analysis.py --configfile /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/Metadata/spsingh_20190521_D1.json --cluster-config peak_analysis.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
#--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
#snakemake --snakefile build_trackhub.py --configfile Encode.json --cluster-config build_trackhub.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
#--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
rm -rf ./*.pos
rm -rf ./*.tmp

echo "Pipeline execution successfully finished at: "$(date)