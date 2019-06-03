Pipeline execution initiated at: Wed May 22 12:34:36 EDT 2019
[+] Loading snakemake  5.4.4 
[+] Loading samtools 1.9  ... 
[+] Loading java 12.0.1  ... 
Building DAG of jobs...
Unlocking working directory.
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 30
Unlimited resources: mem_mb
Job counts:
	count	jobs
	1	End_Point
	16	Pre_Process_Paired
	17

[Wed May 22 12:40:02 2019]
Job 3: Pre_Process_Paired: Ccl20-_UT|AB3825-L1

Submitted job 3 with external jobid '27130916'.

[Wed May 22 12:40:05 2019]
Job 16: Pre_Process_Paired: Ccl20+_T|AB3838-L1

Submitted job 16 with external jobid '27130918'.

[Wed May 22 12:40:07 2019]
Job 6: Pre_Process_Paired: Ccl20+_UT|AB3828-L1

Submitted job 6 with external jobid '27130919'.

[Wed May 22 12:40:10 2019]
Job 12: Pre_Process_Paired: Ccl20-_T|AB3834-L1

Submitted job 12 with external jobid '27130920'.

[Wed May 22 12:40:13 2019]
Job 10: Pre_Process_Paired: Ccl20-_T|AB3832-L1

Submitted job 10 with external jobid '27130921'.

[Wed May 22 12:40:15 2019]
Job 14: Pre_Process_Paired: Ccl20+_T|AB3836-L1

Submitted job 14 with external jobid '27130922'.

[Wed May 22 12:40:18 2019]
Job 11: Pre_Process_Paired: Ccl20-_T|AB3833-L1

Submitted job 11 with external jobid '27130923'.

[Wed May 22 12:40:22 2019]
Job 4: Pre_Process_Paired: Ccl20-_UT|AB3826-L1

Submitted job 4 with external jobid '27130924'.

[Wed May 22 12:40:24 2019]
Job 5: Pre_Process_Paired: Ccl20+_UT|AB3827-L1

Submitted job 5 with external jobid '27130925'.

[Wed May 22 12:40:26 2019]
Job 8: Pre_Process_Paired: Ccl20+_UT|AB3830-L1

Submitted job 8 with external jobid '27130926'.

[Wed May 22 12:40:29 2019]
Job 9: Pre_Process_Paired: Ccl20-_T|AB3831-L1

Submitted job 9 with external jobid '27130927'.

[Wed May 22 12:40:32 2019]
Job 15: Pre_Process_Paired: Ccl20+_T|AB3837-L1

Submitted job 15 with external jobid '27130928'.

[Wed May 22 12:40:35 2019]
Job 1: Pre_Process_Paired: Ccl20-_UT|AB3823-L1

Submitted job 1 with external jobid '27130929'.

[Wed May 22 12:40:37 2019]
Job 13: Pre_Process_Paired: Ccl20+_T|AB3835-L1

Submitted job 13 with external jobid '27130930'.

[Wed May 22 12:40:40 2019]
Job 2: Pre_Process_Paired: Ccl20-_UT|AB3824-L1

Submitted job 2 with external jobid '27130931'.

[Wed May 22 12:40:43 2019]
Job 7: Pre_Process_Paired: Ccl20+_UT|AB3829-L1

Submitted job 7 with external jobid '27130932'.
[Wed May 22 13:04:24 2019]
Finished job 1.
1 of 17 steps (6%) done
[Wed May 22 13:07:08 2019]
Finished job 6.
2 of 17 steps (12%) done
[Wed May 22 13:07:20 2019]
Finished job 2.
3 of 17 steps (18%) done
[Wed May 22 13:08:04 2019]
Finished job 12.
4 of 17 steps (24%) done
[Wed May 22 13:08:26 2019]
Finished job 4.
5 of 17 steps (29%) done
[Wed May 22 13:09:32 2019]
Finished job 3.
6 of 17 steps (35%) done
[Wed May 22 13:14:40 2019]
Finished job 5.
7 of 17 steps (41%) done
[Wed May 22 13:15:11 2019]
Finished job 15.
8 of 17 steps (47%) done
[Wed May 22 13:15:41 2019]
Finished job 7.
9 of 17 steps (53%) done
[Wed May 22 13:17:11 2019]
Finished job 9.
10 of 17 steps (59%) done
[Wed May 22 13:19:01 2019]
Finished job 16.
11 of 17 steps (65%) done
[Wed May 22 13:20:12 2019]
Finished job 13.
12 of 17 steps (71%) done
[Wed May 22 13:21:12 2019]
Finished job 14.
13 of 17 steps (76%) done
[Wed May 22 13:26:42 2019]
Finished job 8.
14 of 17 steps (82%) done
[Wed May 22 13:28:02 2019]
Finished job 11.
15 of 17 steps (88%) done
[Wed May 22 13:31:12 2019]
Finished job 10.
16 of 17 steps (94%) done

[Wed May 22 13:31:12 2019]
localrule End_Point:
    input: /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_UT/pre_process/AB3823-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_UT/pre_process/AB3824-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_UT/pre_process/AB3825-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_UT/pre_process/AB3826-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_UT/pre_process/AB3827-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_UT/pre_process/AB3828-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_UT/pre_process/AB3829-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_UT/pre_process/AB3830-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_T/pre_process/AB3831-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_T/pre_process/AB3832-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_T/pre_process/AB3833-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20-_T/pre_process/AB3834-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_T/pre_process/AB3835-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_T/pre_process/AB3836-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_T/pre_process/AB3837-L1.R1.processed.fastq, /data/shamsaddinisha/jfarber-20160119-P1/spsingh-20190308-E1/spsingh-20190521-D1/mm10/Ccl20+_T/pre_process/AB3838-L1.R1.processed.fastq
    jobid: 0

[Wed May 22 13:31:12 2019]
Finished job 0.
17 of 17 steps (100%) done
Complete log: /gpfs/gsfs6/users/shamsaddinisha/ATAC_Seq/Snakemake/.snakemake/log/2019-05-22T123722.767614.snakemake.log
Pipeline execution successfully finished at: Wed May 22 13:31:22 EDT 2019
