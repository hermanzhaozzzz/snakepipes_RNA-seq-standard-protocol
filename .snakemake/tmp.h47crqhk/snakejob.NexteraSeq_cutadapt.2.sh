#!/bin/bash
# properties = {"type": "single", "rule": "NexteraSeq_cutadapt", "local": false, "input": ["../fastq/pcif1-WT-2_R1.fq.gz", "../fastq/pcif1-WT-2_R2.fq.gz"], "output": ["../fix.fastq/293T-RNASeq-pcif1-WT-2_R1_cutadapt.fq.gz", "../fix.fastq/293T-RNASeq-pcif1-WT-2_R2_cutadapt.fq.gz"], "wildcards": {"sample": "pcif1-WT-2"}, "params": {}, "log": ["../fix.fastq/293T-RNASeq-pcif1-WT-2_cutadapt.log"], "threads": 1, "resources": {}, "jobid": 2, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/2021_Other_projects/2021_CTT/snakepipes_RNA-seq-standard-protocol && \
PATH='/home/zhaohuanan/zhaohn_HD/miniconda3/bin':$PATH /home/zhaohuanan/zhaohn_HD/miniconda3/bin/python3.8 \
-m snakemake ../fix.fastq/293T-RNASeq-pcif1-WT-2_R1_cutadapt.fq.gz --snakefile /gpfs/user/zhaohuanan/3.project/2021_Other_projects/2021_CTT/snakepipes_RNA-seq-standard-protocol/Snakefile.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/user/zhaohuanan/3.project/2021_Other_projects/2021_CTT/snakepipes_RNA-seq-standard-protocol/.snakemake/tmp.h47crqhk ../fastq/pcif1-WT-2_R1.fq.gz ../fastq/pcif1-WT-2_R2.fq.gz --latency-wait 60 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules NexteraSeq_cutadapt --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

