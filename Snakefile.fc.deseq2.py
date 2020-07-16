FEATURECOUNTS = "/home/zhaohuanan/anaconda3/envs/testDE/bin/featureCounts"
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/annotate_hg38/20200714_ComprehensiveGeneAnnotation_Chr_gencode.v29.annotation.gtf"

# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    '2Mock-1_combined',# 前面试bg后面是exp
    '2BE4-All-1_combined',
    '2Vector-1_combined',
#     'BE3-1_combined',
#     'BE3-2_combined',
    'BE4-1_combined',
    'BE4-2_combined',
    'EM-1_combined',
    'EM-2_combined',
    'BE4-0706-rep1',
    'M2-0706-rep1'
]

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        '../featureCounts/all_feature.txt'

rule featureCounts:
#     input:
    params:
        BAM = ' '.join(expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam",sample=SAMPLES))
    output:
        '../featureCounts/all_feature.txt'
    log:
        '../featureCounts/run_FC.log'
    shell:
        """
        {FEATURECOUNTS} \
        -T 24 \
        -p \
        -t exon \
        -g gene_id \
        -a {HG38_GTF} \
        -o {output} \
        {params.BAM} 2>{log}
        """

# -T 使用的线程数
# -p 如果是paird end 就用
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个
