# _*_ coding: UTF-8 _*_

########################################################################
# ZHAO Huanan
# 2020-05-02
# 293T RNASeq analysis pipeline
######################################################################## 
# run on abyss
# /home/zhaohuanan/zhaohn_HD/3.project/2.repeat_David_Liu_2020_NBT/3.test_snake/raw_bam

# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# before this, makesure you done fastqc + multiqc and know how to cutadapt and trim
# 1. cutadapt
# 2. STAR alingment 
# 3. picard change RG
# 4. calculate FPKM
# 5. samtools sort by position
# 6. picard mark duplicates
# 7. build index

# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
# cutadapt --v 1.18
CUTADAPT = "/home/zhaohuanan/zhaohn_HD/1.apps/anaconda3/envs/py37/bin/cutadapt"
# STAR --version 2.7.2b
STAR = "/gpfs/build/bin/STAR"
# JAVA version 1.8
# openjdk 11.0.1 2018-10-16 LTS
# OpenJDK Runtime Environment Zulu11.2+3 (build 11.0.1+13-LTS)
# OpenJDK 64-Bit Server VM Zulu11.2+3 (build 11.0.1+13-LTS, mixed mode)
JAVA = "/home/zhaohuanan/zhaohn_HD/1.apps/anaconda3/envs/py37/bin/java"
PICARD = "/gpfs/user/zhaohuanan/1.apps/picard/picard.jar"
# cufflinks v2.2.1
CUFFLINKS = "/home/zhaohuanan/zhaohn_HD/1.apps/anaconda3/envs/py27/bin/cufflinks"
# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
HG38_FA_DICT = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
# STAR --version 2.7.2b
STAR_HG38_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/star_hg38"
HG39_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/annotate_hg38/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf"
# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
#      "BE3-1_combined",
#      "BE3-2_combined",
#      "BE4-1_combined",
#      "BE4-2_combined",
#      "EM-1_combined",
#      "EM-2_combined"
    "test_combined",
    "test2_combined"
]
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",sample=SAMPLES),
        expand("../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned.out.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam",sample=SAMPLES),
        expand("../fpkm/{sample}",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned.out.bam.bai",sample=SAMPLES)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cutadapter
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule TruSeq_cutadapt:
    input:
        "../reads/{sample}_R1.fastq.gz",
        "../reads/{sample}_R2.fastq.gz"
    output:
        "../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",
        "../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz"
    log:
        "../fix.fastq/293T-RNASeq-{sample}_cutadapt.log"
    shell:# using illumina universal adaptor
        "srun -T 24 -c 24 \
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 55 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1 "    

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# STAR mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule STAR_mapping:
    input:
        fq1 = "../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",
        fq2 = "../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    log:
        "../bam/293T-RNASeq-{sample}_Aligned.out.log"
    params:
        "../bam/293T-RNASeq-{sample}_"
    shell:
        "srun -T 24 {STAR} \
        --genomeDir {STAR_HG38_INDEX} \
        --runThreadN 24 \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params} \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical > {log} 2>&1"
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# add @RG tag (mostly for GATK SNP/SNV calling)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule add_RG_tag:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam"
    params:
        tag = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        "srun -T 24 samtools addreplacerg -r {params.tag} -@ 24 -O BAM -o {output} --reference {HG38_FA_DICT} {input}"
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cufflinks calculate FPKM
# ------------------------------------------------------------------------------------------>>>>>>>>>>     
rule cufflinks_FPKM:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    output:
        directory('../fpkm/{sample}')
    shell:
        """
        srun -T 24 -c 24 \
        {CUFFLINKS} -p 24 --library-type fr-firststrand \
        -G {HG39_GTF} \
        -o {output} \
        {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam"
    shell:
        "srun -T 24 -c 24 samtools sort -O BAM -o {output} -T {output}.temp -@ 24 -m 4G {input}"    
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# picard mark duplicate
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_mark_duplicate:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.bam",
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.matrix"
    log:
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.log"
    shell:
        "srun -T 24 -c 24 {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 -jar {PICARD} MarkDuplicates I={input} O={output[0]} M={output[1]} ASO=coordinate 2>{log}"

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>      
rule BAM_index:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam.bai"
    shell:
        "srun -T 24 samtools index -@ 24 {input} {output}"   





# ##############################
# # SplitNCigarReads for RNA(DNA no need)
# ##############################
# rule GATK_SplitNCigarReads:
#     input:
#         "bam/{sample}_Aligned_sortn_rmDup.bam"
#     output:
#         "bam/{sample}_Aligned_sortn_SplitNCigar.bam"
#     log:
#         "bam/{sample}_Aligned_sortn_SplitNCigar.bam.log"
#     shell:
#         """
#         srun -T 24 -c 24 \
#         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
#         -jar {GATK4} \
#         SplitNCigarReads \
#         -I {input} \
#         -O {output} \
#         -R {HG38_FA} 2>{log}
#         """ 
        
        
        
        
    

# # ##############################
# # # HaplotypeCaller -> g.vcf
# # ##############################
# # rule GATK_HaplotypeCaller:
# #     input:
# #         "bam/{sample}_Aligned_sort_MarkDup.bam"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_raw_variants.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_raw_variants.vcf.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         HaplotypeCaller \
# #         -I {input} \
# #         -O {output} \
# #         -R {HG38_FA} 2>{log}
# #         """ 
# # ##############################
# # # SelectVariants SNP
# # ############################## 
# # rule GATK_SelectVariants_SNP:
# #     input:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_raw_variants.vcf"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_SNP.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_SNP.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         SelectVariants \
# #         -select-type SNP \
# #         -V {input} \
# #         -O {output} 2>{log}
# #         """

# # ##############################
# # # VariantFiltration SNP
# # ##############################
# # rule GATK_VariantFiltration_SNP:
# #     input:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_SNP.vcf"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_SNP.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_SNP.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         VariantFiltration \
# #         --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
# #         --filter-name "PASS" \
# #         -V {input} \
# #         -O {output} 2>{log}
# #         """
# # ##############################
# # # SelectVariants INDEL
# # ##############################
# # rule GATK_SelectVariants_INDEL:
# #     input:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_raw_variants.vcf"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_INDEL.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_INDEL.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         SelectVariants \
# #         -select-type INDEL \
# #         -V {input} \
# #         -O {output} 2>{log}
# #         """
# # ##############################
# # # VariantFiltration INDEL
# # ##############################
# # rule GATK_VariantFiltration_INDEL:
# #     input:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_select_variants_INDEL.vcf"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_INDEL.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_INDEL.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         VariantFiltration \
# #         --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
# #         --filter-name "PASS" \
# #         -V {input} \
# #         -O {output} 2>{log}
# #         """
# # ##############################
# # # MergeVcfs -> vcf
# # ##############################
# # rule GATK_MergeVcfs:
# #     input:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_SNP.vcf",
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_INDEL.vcf"
# #     output:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_merge.vcf"
# #     log:
# #         "vcf/{sample}_gatk4_HaplotypeCaller_filtered_variants_merge.log"
# #     shell:
# #         """
# #         srun -T 24 -c 24 \
# #         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads=24 \
# #         -jar {GATK4} \
# #         MergeVcfs \
# #         -I {input[0]} \
# #         -I {input[1]} \
# #         -O {output} 2>{log}
# #         """