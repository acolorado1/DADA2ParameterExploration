rule all:
    input:
        "output/TrimInfo.tsv",
        "output/sumEEInfo.tsv",
        "output/ScatterReadLengthByAvgEE.png",
        "output/HeatmapIndexValueByAvgEE.png",
        "output/HistogramRetainedReadCount.png"

rule create_TSVs:
    input:
        "LOZ-CSU-Nano_S1_L001_R1_001.fastq"
    output:
        TrimInfo = "output/TrimInfo.tsv",
        SumEEInfo = "output/sumEEInfo.tsv"
    params:
        threshold = 30,
        o_mtp = 0.1,
        a_mtp = 0.2
    shell:
        "python " + "smartdada2/main.py --fq {input} --t_of {output.TrimInfo} --th {params.threshold} --o_mtp {params.o_mtp} --a_mtp {params.a_mtp} --EE_of {output.SumEEInfo}"

rule get_plots:
    input:
        TrimInfo = "output/TrimInfo.tsv",
        sumEEInfo = "output/sumEEInfo.tsv"
    output:
        "output/ScatterReadLengthByAvgEE.png",
        "output/HeatmapIndexValueByAvgEE.png",
        "output/HistogramRetainedReadCount.png"
    shell:
        "Rscript " + "smartdada2/R_scripts/vizualize.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo}"
