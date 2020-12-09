rule STAR_notrim:
    input:
        fwd = "samples/raw/{sample}_R1.fastq.gz",
        rev = "samples/raw/{sample}_R2.fastq.gz"
    output:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star_notrim/{sample}_bam/Log.final.out",
        "samples/star_notrim/{sample}_bam/ReadsPerGene.out.tab"
    threads: 12
    params:
        gtf=config["gtf_file"]
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star_notrim/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --twopassMode Basic
                """)


rule samtools_index:
    input:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
    output:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam.bai",
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        "samtools index {input} {output}"


rule star_statistics:
    input:
        expand("samples/star_notrim/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"


rule samtools_stats:
    input:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/samtools_stats/{sample}.txt"
    conda:
        "../envs/omic_qc_wf.yaml"
    wrapper:
        "0.17.0/bio/samtools/stats"


rule compile_star_counts:
    input:
        expand("samples/star_notrim/{sample}_bam/ReadsPerGene.out.tab",sample=SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_star_counts.py"

