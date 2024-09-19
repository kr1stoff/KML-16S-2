rule plot_otu_venn:
    input:
        otu_tsv=rules.diversity_core_export.output,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        img_dir=directory('diff/Venn/otu')
    benchmark:
        '.log/diff/plot_otu_venn.bm'
    log:
        '.log/diff/plot_otu_venn.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output.img_dir}
        Rscript {config[my_scripts]}/diff_otu_venn.R {input.otu_tsv} {input.metadata} group {output.img_dir}/group.otu.venn > {log} 2>&1
        """


rule plot_taxa_venn:
    input:
        table_tsv=rules.taxa_collapse.output.table_tsv,
        metadata=rules.select_sampling_depth.output.re_metadata,
    output:
        directory('diff/Venn/{levels}')
    benchmark:
        '.log/diff/plot_taxa_venn_{levels}.bm'
    log:
        '.log/diff/plot_taxa_venn_{levels}.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_otu_venn.R {input.table_tsv} {input.metadata} group {output}/group.{wildcards.levels}.venn > {log} 2>&1
        """


rule plot_otu_PCA:
    input:
        otu_tsv=rules.diversity_core_export.output,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/PCA/otu')
    benchmark:
        '.log/diff/plot_otu_PCA.bm'
    log:
        '.log/diff/plot_otu_PCA.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_otu_PCA.R {input.otu_tsv} {input.metadata} group {output}/group.otu.PCA "OTU levels" > {log} 2>&1
        """


rule plot_taxa_PCA:
    input:
        table_tsv=rules.taxa_collapse.output.table_tsv,
        metadata=rules.select_sampling_depth.output.re_metadata,
    output:
        directory('diff/PCA/{levels}'),
    benchmark:
        '.log/diff/plot_taxa_PCA_{levels}.bm'
    log:
        '.log/diff/plot_taxa_PCA_{levels}.log'
    conda:
        config['conda']['microplot']
    shell:
        """        
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_otu_PCA.R {input.table_tsv} {input.metadata} group \
            {output}/group.{wildcards.levels}.PCA "{wildcards.levels} levels" > {log} 2>&1
        """


rule plot_otu_PLSDA:
    input:
        otu_tsv=rules.diversity_core_export.output,
        metadata=rules.select_sampling_depth.output.re_metadata,
    output:
        directory('diff/PLS-DA/otu')
    benchmark:
        '.log/diff/plot_otu_PLSDA.bm'
    log:
        '.log/diff/plot_otu_PLSDA.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_PLS-DA.R {input.otu_tsv} {input.metadata} group {output}/group.otu.PLS-DA "OTU" > {log} 2>&1
        """


rule plot_taxa_PLSDA:
    input:
        table_tsv=rules.taxa_collapse.output.table_tsv,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/PLS-DA/{levels}')
    benchmark:
        '.log/diff/plot_taxa_PLSDA_{levels}.bm'
    log:
        '.log/diff/plot_taxa_PLSDA_{levels}.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_PLS-DA.R {input.table_tsv} {input.metadata} group \
            {output}/group.{wildcards.levels}.PLS-DA "{wildcards.levels}" > {log} 2>&1
        """


rule plot_otu_OPLSDA:
    input:
        otu_tsv=rules.diversity_core_export.output,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/OPLS-DA/otu')
    benchmark:
        '.log/diff/plot_otu_OPLSDA.bm'
    log:
        '.log/diff/plot_otu_OPLSDA.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_OPLS-DA.R {input.otu_tsv} {input.metadata} group {output}/group.otu.OPLS-DA "OTU" > {log} 2>&1
        # 因为拼图的包，可能在工作目录生成 Rplots.pdf
        if [ -f Rplots.pdf ];then rm Rplots.pdf;fi
        """


rule plot_taxa_OPLSDA:
    input:
        table_tsv=rules.taxa_collapse.output.table_tsv,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/OPLS-DA/{levels}')
    benchmark:
        '.log/diff/plot_taxa_OPLSDA_{levels}.bm'
    log:
        '.log/diff/plot_taxa_OPLSDA_{levels}.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_OPLS-DA.R {input.table_tsv} {input.metadata} group \
            {output}/group.{wildcards.levels}.OPLS-DA "{wildcards.levels}" > {log} 2>&1
        if [ -f Rplots.pdf ];then rm Rplots.pdf;fi
        """


rule plot_RandomForest:
    input:
        otu_tsv=rules.diversity_core_export.output,
        taxonomy_dir=rules.taxa.output.taxonomy_dir,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/RandomForest')
    benchmark:
        '.log/diff/plot_RandomForest.bm'
    log:
        '.log/diff/plot_RandomForest.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_random_forest.R {input.otu_tsv} {input.taxonomy_dir}/taxonomy.tsv \
            {input.metadata} group {output}/group.RandomForest > {log} 2>&1
        """


rule plot_Lefse:
    input:
        otu_tsv=rules.diversity_core_export.output,
        taxonomy_dir=rules.taxa.output.taxonomy_dir,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('diff/Lefse')
    benchmark:
        '.log/diff/plot_Lefse.bm'
    log:
        '.log/diff/plot_Lefse.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/diff_lefse.R {input.otu_tsv} {input.taxonomy_dir}/taxonomy.tsv {input.metadata} \
                group {output}/group.Lefse > {log} 2>&1
        """
