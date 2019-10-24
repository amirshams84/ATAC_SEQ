#!/bin/sh
# properties = {"type": "single", "rule": "Pre_Process_Paired", "local": false, "input": [], "output": ["/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Encode/pre_process/ctl1.R1.processed.fastq", "/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Encode/pre_process/ctl1.R2.processed.fastq"], "wildcards": {"design": "Encode", "sample": "ctl1"}, "params": {}, "log": [], "threads": 10, "resources": {"mem_mb": 10000}, "jobid": 1, "cluster": {"partition": "quick", "time": 240, "jobname": "Snake.Pre_Process_Paired.design=Encode,sample=ctl1", "extra": "--gres=lscratch:100", "output": "./logs/pre_process/Pre_Process_Paired.design=Encode,sample=ctl1.stdout", "error": "./logs/pre_process/Pre_Process_Paired.design=Encode,sample=ctl1.stderr"}}
cd /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake && \
/usr/local/Anaconda/envs_app/snakemake/5.4.4/bin/python3.6 \
-m snakemake '/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Encode/pre_process/ctl1.R1.processed.fastq' --snakefile /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/pre_process.py \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.ckzn9uiq --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/Metadata/spsingh_20190513_D4.json  --allowed-rules Pre_Process_Paired --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.ckzn9uiq/1.jobfinished" || (touch "/gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.ckzn9uiq/1.jobfailed"; exit 1)

