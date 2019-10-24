#!/bin/sh
# properties = {"type": "single", "rule": "Pre_Process_Paired", "local": false, "input": ["/data/RTB/GTS/Solexa/Bcl2fastq/190510_NS500189_0127_AHMWKVBGX7/spsingh-20190425-H1/AB3826-L1_HMWKVBGX7_S4_R1_001.fastq.gz"], "output": ["/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Ccl20-_UT4/pre_process/AB3826-L1.R1.processed.fastq", "/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Ccl20-_UT4/pre_process/AB3826-L1.R2.processed.fastq"], "wildcards": {"design": "Ccl20-_UT4", "sample": "AB3826-L1"}, "params": {}, "log": [], "threads": 20, "resources": {"mem_mb": 100000}, "jobid": 4, "cluster": {"partition": "norm", "time": 720, "jobname": "Snake.Pre_Process_Paired.design=Ccl20-_UT4,sample=AB3826-L1", "extra": "--gres=lscratch:100", "output": "./logs/pre_process/Pre_Process_Paired.design=Ccl20-_UT4,sample=AB3826-L1.stdout", "error": "./logs/pre_process/Pre_Process_Paired.design=Ccl20-_UT4,sample=AB3826-L1.stderr"}}
cd /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake && \
/usr/local/Anaconda/envs_app/snakemake/5.4.4/bin/python3.6 \
-m snakemake '/data/shamsaddinisha/ jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190513-D4/mm10/Ccl20-_UT4/pre_process/AB3826-L1.R1.processed.fastq' --snakefile /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/pre_process.py \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.i7_dlm6c /data/RTB/GTS/Solexa/Bcl2fastq/190510_NS500189_0127_AHMWKVBGX7/spsingh-20190425-H1/AB3826-L1_HMWKVBGX7_S4_R1_001.fastq.gz --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/Metadata/spsingh_20190513_D4.json  --allowed-rules Pre_Process_Paired --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.i7_dlm6c/4.jobfinished" || (touch "/gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/tmp.i7_dlm6c/4.jobfailed"; exit 1)

