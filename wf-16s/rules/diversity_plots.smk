rule plot_rare_curve:
    input:
        rules.alpha_rarefaction_export.output
    output:
        directory('diversity/alpha_rare_curve')
    benchmark:
        '.log/diversity/plot_rare_curve.bm'
    log:
        '.log/diversity/plot_rare_curve.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/alpha_rare_curve.R {input}/observed_features.csv group {output}/alpha_rare_cure> {log} 2>&1
        """


rule plot_specaccum:
    input:
        rules.dada2_rename_table.output.re_tsv
    output:
        'diversity/alpha_specaccum/specaccum.{fmt}'
    benchmark:
        '.log/diversity/plot_specaccum_{fmt}.bm'
    log:
        '.log/diversity/plot_specaccum_{fmt}.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        Rscript {config[my_scripts]}/alpha_specaccum.R {input} {output} > {log} 2>&1
        """


rule plot_rank_abundance_curve:
    input:
        rules.diversity_core_export.output
    output:
        directory('diversity/alpha_rank_abundance_curve')
    benchmark:
        '.log/diversity/plot_rank_abundance_curve.bm'
    log:
        '.log/diversity/plot_rank_abundance_curve.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/alpha_rank_abundance_curve.R {input} {output}/alpha_rank_abundance > {log} 2>&1
        """


rule plot_beta_PCOA:
    input:
        asv_tsv=rules.diversity_core_export.output,
        species_tsv='taxa/collapse/Species/feature-table.tsv',
        genus_tsv='taxa/collapse/Genus/feature-table.tsv',
        family_tsv='taxa/collapse/Family/feature-table.tsv',
        core_dir=rules.diversity_core.output.metrics_dir
    output:
        directory('diversity/beta_PCOA')
    benchmark:
        '.log/diversity/plot_beta_PCOA.bm'
    log:
        '.log/diversity/plot_beta_PCOA.log'
    conda:
        config['conda']['microplot']
    params:
        metadata=config['metadata']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} {input.asv_tsv} \
            {input.core_dir}/bray_curtis_distance_matrix/distance-matrix.tsv group bray {output}/beta_PCOA.ASV.bray \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} {input.asv_tsv} \
            {input.core_dir}/jaccard_distance_matrix/distance-matrix.tsv group jaccard {output}/beta_PCOA.ASV.jaccard \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} {input.asv_tsv} \
            {input.core_dir}/weighted_unifrac_distance_matrix/distance-matrix.tsv group euclidean {output}/beta_PCOA.ASV.weighted_unifrac \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} {input.asv_tsv} \
            {input.core_dir}/unweighted_unifrac_distance_matrix/distance-matrix.tsv group manhattan {output}/beta_PCOA.ASV.unweighted_unifrac \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} \
            taxa/collapse/Species/feature-table.tsv none group bray {output}/beta_PCOA.Species.bray \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} \
            taxa/collapse/Genus/feature-table.tsv none group bray {output}/beta_PCOA.Genus.bray \
            2>> {log}
        Rscript {config[my_scripts]}/beta_PCOA_by_matrix.R {params.metadata} \
            taxa/collapse/Family/feature-table.tsv none group bray {output}/beta_PCOA.Family.bray \
            2>> {log}
        """


rule plot_beta_NMDS:
    input:
        rules.select_sampling_depth.output.re_metadata,
        rules.diversity_core.output.metrics_dir,
        # 需要确保该 rule 先完成
        rules.diversity_core_export.output
    output:
        directory('diversity/beta_NMDS')
    benchmark:
        '.log/diversity/plot_beta_NMDS.bm'
    log:
        '.log/diversity/plot_beta_NMDS.log'
    threads:
        config['threads']['low']
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        for method in bray_curtis jaccard unweighted_unifrac weighted_unifrac
        do
            echo "Rscript {config[my_scripts]}/beta_NMDS.R {input[0]} \
                   {input[1]}/${{method}}_distance_matrix/distance-matrix.tsv group {output}/${{method}}.group"
        done | parallel -j {threads} > {log} 2>&1
        """


rule plot_beta_UPGMA:
    input:
        rules.select_sampling_depth.output.re_metadata,
        rules.diversity_core.output.metrics_dir,
        rules.diversity_core_export.output
    output:
        directory('diversity/beta_UPGMA')
    benchmark:
        '.log/diversity/plot_beta_UPGMA.bm'
    log:
        '.log/diversity/plot_beta_UPGMA.log'
    threads:
        config['threads']['low']
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        for method in bray_curtis jaccard unweighted_unifrac weighted_unifrac
        do
            echo "Rscript {config[my_scripts]}/beta_UPGMA.R {input[1]}/${{method}}_distance_matrix/distance-matrix.tsv \
                {input[0]} group {output}/${{method}}.group"
        done | parallel -j {threads} > {log} 2>&1
        """


rule plot_beta_heatmap:
    input:
        rules.select_sampling_depth.output.re_metadata,
        rules.diversity_core.output.metrics_dir,
        rules.diversity_core_export.output
    output:
        directory('diversity/beta_heatmap')
    benchmark:
        '.log/diversity/plot_beta_heatmap.bm'
    log:
        '.log/diversity/plot_beta_heatmap.log'
    threads:
        config['threads']['low']
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        for method in bray_curtis jaccard unweighted_unifrac weighted_unifrac
        do
            echo "Rscript {config[my_scripts]}/beta_heatmap.R {input[1]}/${{method}}_distance_matrix/distance-matrix.tsv \
                {input[0]} group {output}/${{method}}.group.pdf"
            echo "Rscript {config[my_scripts]}/beta_heatmap.R {input[1]}/${{method}}_distance_matrix/distance-matrix.tsv \
                {input[0]} group {output}/${{method}}.group.png"
        done | parallel -j {threads} > {log} 2>&1
        """


rule beta_Adonis:
    input:
        otu_tsv=rules.diversity_core_export.output,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diversity/beta_Adonis')
    benchmark:
        '.log/diversity/beta_Adonis.bm'
    log:
        '.log/diversity/beta_Adonis.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/beta_Adonis.R {input.otu_tsv} {input.metadata} {output}/beta_Adonis.group.tmp group > {log} 2>&1

        #根据首列去重, 实际只是去除了重复的表头
        cat {output}/beta_Adonis.*.tmp | \
            awk '!a[$1]++' > {output}/beta_Adonis.tsv && \
            rm -rf {output}/beta_Adonis.*.tmp >> {log} 2>&1
        """
