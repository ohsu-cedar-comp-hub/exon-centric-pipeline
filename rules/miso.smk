rule generate_map:
    params:
        events = config['gff'],
    output:
        "samples/miso/settings/{type}_map.pkl"
    script:
        "../scripts/generate_map_file.py"


rule compute_insert_length:
    input:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
    params:
        const_exons = config['const_exons']
    output:
        "samples/miso/analysis/{sample}/insert-dist/Aligned.sortedByCoord.out.bam.insert_len"
    conda:
        "../envs/miso_bioconda.yaml"
    shell:
        """
        pe_utils --compute-insert-len {input} {params.const_exons} --output-dir samples/miso/analysis/{wildcards.sample}/insert-dist
        """


rule run_miso:
    input:
        "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/miso/analysis/{sample}/insert-dist/Aligned.sortedByCoord.out.bam.insert_len"
    params:
        events = config['events'],
        read_len = config['read_length']
    conda:
        "../envs/miso_bioconda.yaml"
    output:
        "samples/miso/analysis/{sample}/{sample}/logs/{sample}_run.txt"
    script:
        "../scripts/run_miso.py"


rule summarize_miso:
    input:
        "samples/miso/analysis/{sample}/{sample}/logs/{sample}_run.txt"
    output:
        "samples/miso/analysis/{sample}/summary/{sample}.miso_summary"
    conda:
        "../envs/miso_bioconda.yaml"
    script:
        "../scripts/summarize_miso.py"


rule compare_miso:
  input:
      sample1 = "samples/miso/analysis/{sample1}/{sample1}",
      sample2 = "samples/miso/analysis/{sample2}/{sample2}"
  output:
      "samples/miso/comparisons/{sample1}_vs_{sample2}/bayes-factors/{sample1}_vs_{sample2}.miso_bf"
  conda:
      "../envs/miso_bioconda.yaml"
  shell:
    """
    compare_miso --compare-samples {input.sample1} {input.sample2} samples/miso/comparisons
    """

rule summarizeVs:
    input:
        sample_list=expand("samples/miso/comparisons/{sample1}_vs_{sample2}",sample1=SAMPLES,sample2=SAMPLES[::-1]),
        map_file="samples/miso/settings/{type}_map.pkl".format(type=config['event_type'])
    params:
        type=config['event_type']
    output:
        expand("samples/miso/summary_vs/{sample1}_vs_{sample2}.miso_bf",sample1=SAMPLES,sample2=SAMPLES[::-1]),
    script:
        "../scripts/summarize_vs.py"


rule link_summaries:
    input:
        "samples/miso/analysis/{sample}/summary/{sample}.miso_summary"
    output:
        "samples/miso/summary_vs/{sample}.miso_summary"
    shell:
        """
        ln -s $(pwd)/{input} {output} 
        """


rule link_bams:
    input:
        bam = "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
        index = "samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam.bai",
    output:
        bam_out = "samples/miso/bams/{sample}",
        index_out = "samples/miso/bams/{sample}.bai"
    shell:
        """
        ln -s $(pwd)/{input.bam} {output.bam_out}
        ln -s $(pwd)/{input.index} {output.index_out}
        """


rule format_sashimi_settings:
    input:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"]),
    params:
        config['omic_meta_data'],
        contrast=get_contrast,
        condition = config['pca']['labels']
    output:
        "samples/miso/settings/sashimi_plot_{contrast}.txt",
    script:
        "../scripts/format_sashimi_settings.py"


rule generate_consolidated_summaries:
    input:
        comparisons=expand("samples/miso/summary_vs/{sample1}_vs_{sample2}.miso_bf",sample1=SAMPLES, sample2=SAMPLES[::-1]),
        summaries=expand("samples/miso/summary_vs/{sample}.miso_summary",sample=SAMPLES),
    params:
        samples=expand("{sample}",sample=SAMPLES),
    output:
        "samples/miso/results/consolidated_summaries.txt"
    script:
        "../scripts/consolidate_summaries.py"


rule generate_group_files:
    params:
        contrast=get_contrast,
        meta = config['omic_meta_data'],
        condition = config['pca']['labels']
    output:
        "samples/miso/group_files/{contrast}_group.txt"
    script:
        "../scripts/make_groups.py"


rule generate_monotonic_file:
    input:
        summary = "samples/miso/results/consolidated_summaries.txt",
        group_f = "samples/miso/group_files/{contrast}_group.txt"
    output:
        "samples/miso/results/monotonic_summaries_{contrast}.txt",
    params:
        config['omic_meta_data'],
        contrast = get_contrast,
        minbayes = config['minbayes'],
        shuffles = config['shuffles']
    conda:
        "../envs/miso_bioconda.yaml"
    script:
        "../scripts/run_monotonic.py"


rule plot_sashimi:
    input:
        mono_fh = "samples/miso/results/monotonic_summaries_{contrast}.txt",
        ps_fh = "samples/miso/settings/sashimi_plot_{contrast}.txt",
        bai = expand("samples/miso/bams/{sample}.bai",sample=SAMPLES)
    output:
        "samples/miso/sashimi_plots/{contrast}/all_plots_all_launched.txt"
    params:
        events = config['events'],
        contrast = get_contrast
    conda:
        "../envs/miso_bioconda.yaml"
    script:
        "../scripts/plot_sashimi_results.py"


rule generate_psi_table:
    input:
        expand("samples/miso/summary_vs/{sample1}_vs_{sample2}.miso_bf",sample1=SAMPLES,sample2=SAMPLES[::-1]),
    output:
        "samples/miso/results/{project_id}_psi_table.txt".format(project_id=config['project_id'])
    script:
        "../scripts/generate_psi_table.py"

