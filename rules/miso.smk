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
        set +u
        source deactivate
        source activate miso_bioconda
        pe_utils --compute-insert-len {input} {params.const_exons} --output-dir samples/miso/analysis/{wildcards.sample}/insert-dist
        """


rule run_miso:
    input:
        bam="samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam",
        insert_len="samples/miso/analysis/{sample}/insert-dist/Aligned.sortedByCoord.out.bam.insert_len",
        index="samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    params:
        events = config['events'],
        read_len = config['read_length'],
        out_dir = lambda w: "samples/miso/analysis/{sample}/{sample}".format(sample=w.sample)
    output:
        "samples/miso/analysis/{sample}/{sample}/logs/{sample}_run.txt",
    shell:
        """
        set +u
        source deactivate
        source activate miso_bioconda
        MEAN=$(cat Aligned.sortedByCoord.out.bam.insert_len | head -n 1 | grep -o -E 'mean=[0-9.]+'|tr 'mean=' '\n')
        SDEV=$(cat Aligned.sortedByCoord.out.bam.insert_len | head -n 1 | grep -o -E 'sdev=[0-9.]+'|tr 'sdev=' '\n')
        miso --run {params.events} {input.bam} --output-dir {params.out_dir} --read-len {params.read_len} --paired-end $MEAN $SDEV
        touch {output}
        """

rule summarize_miso:
    input:
        "samples/miso/analysis/{sample}/{sample}/logs/{sample}_run.txt"
    output:
        "samples/miso/analysis/{sample}/summary/{sample}.miso_summary"
    params:
        out_dir = lambda w: "samples/miso/analysis/{sample}".format(sample=w.sample),
        in_dir = lambda w: "samples/miso/analysis/{sample}/{sample}".format(sample=w.sample)
    shell:
        """
        set +u
        source deactivate
        source activate miso_bioconda
        summarize_miso --summarize-samples {params.in_dir} {params.out_dir}
        """


rule compare_miso:
  input:
    sample1 ="samples/miso/analysis/{sample1}/{sample1}/logs/{sample1}_run.txt",
    sample2 = "samples/miso/analysis/{sample2}/{sample2}/logs/{sample2}_run.txt"
  output:
      "samples/miso/comparisons/{sample1}_vs_{sample2}/bayes-factors/{sample1}_vs_{sample2}.miso_bf",
  params:
        smp1 = lambda w: "samples/miso/analysis/{}/{}".format(w.sample1,w.sample1),
        smp2 = lambda w: "samples/miso/analysis/{}/{}".format(w.sample2,w.sample2)
  shell:
    """
    set +u
    source deactivate
    source activate miso_bioconda
    compare_miso --compare-samples {params.smp1} {params.smp2} samples/miso/comparisons
    """

rule summarizeVs:
    input:
        sample_list="samples/miso/comparisons/{sample1}_vs_{sample2}/bayes-factors/{sample1}_vs_{sample2}.miso_bf",
        map_file="samples/miso/settings/{type}_map.pkl".format(type=config['event_type'])
    params:
        type=config['event_type'],
    output:
        "samples/miso/summary_vs/{sample1}_vs_{sample2}.miso_bf"
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
        comparisons=expand("samples/miso/summary_vs/{sample1}_vs_{sample2}.miso_bf",sample1=comb1, sample2=comb2),
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

