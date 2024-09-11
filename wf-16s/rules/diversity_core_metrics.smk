rule select_sampling_depth:
    input:
        table_tsv=rules.dada2_rename_table.output.re_tsv,
        stats=rules.stats_export.output
    output:
        sample_depth_info='diversity/depth/sample_depth_info.json',
        selected_depth='diversity/depth/selected.depth',
        re_metadata='diversity/depth/re-sample-metadata.tsv'
    benchmark:
        '.log/diversity/select_sampling_depth.bm'
    log:
        '.log/diversity/select_sampling_depth.log'
    conda:
        config['conda']['qiime2']
    params:
        '-p 0.9 -d 0.0001'
    shell:
        """
        python {config[my_scripts]}/select_sampling_depth_v2.py \
            {params} \
            -f {input.table_tsv} \
            -o {output.sample_depth_info} > {log} 2>&1

        depth=$(grep "p_sampling_depth" {output.sample_depth_info} | awk '{{print $2}}' | cut -d "," -f 1)
        echo -e "depth:\\t$depth" > {output.selected_depth}
        
        # 生成删除小于抽平取样 depth 的 metadata
        python {config[my_scripts]}/filtered_metadata_samples.py \
            {config[metadata]} \
            {input.stats}/stats.tsv \
            $depth > {output.re_metadata} 2>> {log}
        """


rule diversity_core:
    input:
        tree=rules.phylogeny.output.rooted_tree,
        table=rules.dada2_rename_table.output.re_table,
        selected_depth=rules.select_sampling_depth.output.selected_depth
    output:
        metrics_dir=directory('diversity/core-metrics'),
        rarefied_table='diversity/rarefied_table.qza'
    benchmark:
        '.log/diversity/diversity_core.bm'
    log:
        '.log/diversity/diversity_core.log'
    conda:
        config['conda']['qiime2']
    threads:
        config['threads']['high']
    params:
        '--m-metadata-file ' + config['metadata']
    shell:
        """
        depth=$(cat {input.selected_depth}| cut -f 2)

        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.tree} \
            --i-table {input.table} \
            --p-sampling-depth $depth \
            --p-n-jobs-or-threads {threads} \
            {params} \
            --output-dir {output.metrics_dir} > {log} 2>&1

        #把抽平后的特征表拷贝到另外单独路径, 方便后续rules调用
        cp {output.metrics_dir}/rarefied_table.qza {output.rarefied_table}
        """


rule diversity_core_export:
    input:
        rules.diversity_core.output.metrics_dir
    output:
        rarefied_table_tsv='feature/rarefied_table/feature-table.tsv',
        doneflag=touch('diversity/core-metrics/.export.done')
    benchmark:
        '.log/diversity/diversity_core_export.bm'
    log:
        '.log/diversity/diversity_core_export.log'
    threads:
        config['threads']['high']
    conda:
        config['conda']['qiime2']
    shell:
        """
        for artifact in $(ls {input}/*qza {input}/*.qzv)
        do
            echo "qiime tools export  --input-path $artifact --output-path ${{artifact%.*}}"
        done | parallel -j {threads} > {log} 2>&1

        #把抽平后的rarefied_table/feature-table.biom 转换为tsv
        biom convert -i {input}/rarefied_table/feature-table.biom \
            -o {output.rarefied_table_tsv} --to-tsv >> {log} 2>&1
        """


rule diversity_rarefied_table_freq:
    input:
        rules.diversity_core.output.rarefied_table
    output:
        rarefied_table_freq='diversity/rarefied_table_freq.qza',
        rarefied_table_freq_dir=directory('diversity/rarefied_table_freq_dir'),
        rarefied_table_freq_tsv='diversity/rarefied_table_freq/feature-table.tsv'
    benchmark:
        '.log/diversity/diversity_core_export.bm'
    log:
        '.log/diversity/diversity_core_export.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime feature-table relative-frequency \
            --i-table {input} \
            --o-relative-frequency-table {output.rarefied_table_freq} > {log} 2>&1 && \
        qiime tools export \
            --input-path {output.rarefied_table_freq} \
            --output-path {output.rarefied_table_freq_dir} >> {log} 2>&1 && \
        biom convert \
            -i {output.rarefied_table_freq_dir}/feature-table.biom \
            -o {output.rarefied_table_freq_tsv} --to-tsv >> {log} 2>&1
        """
