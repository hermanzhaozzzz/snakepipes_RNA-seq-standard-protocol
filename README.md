1. https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
先跑质控流程，主要是看看需不需要trim 5‘ 端的不稳定序列和3'端质量不够高的reads
2. 跑mapping流程，如果是star mapping，就用这个流程即可
3. 接下来可以做如GATK SNP calling和差异表达分析