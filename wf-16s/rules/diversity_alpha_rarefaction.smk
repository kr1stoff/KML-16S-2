rule alpha_rarefaction:
    input:
        tree=rules.phylogeny.output.rooted_tree,
        table=rules.dada2_rename_table.output.re_table,
        stats=rules.stats_export.output
    output:
        'diversity/alpha-rarefaction.qzv'
    benchmark:
        '.log/diversity/alpha_rarefaction.bm'
    log:
        '.log/diversity/alpha_rarefaction.log'
    conda:
        config['conda']['qiime2']
    params:
        metadata='--m-metadata-file ' + config['metadata'],
        p_steps=25,
        p_iters=10
    shell:
        """
        echo -n > {log}
        # 获取样本reads值的中位数值, 注意stats前2行是表头和注释, 第6列是过滤后的reads
        # 自动获取样本数据量的中位数值作为p_max_depth
        TargetDepth=$(sort -n -k 6 {input.stats}/stats.tsv | \
            awk -v total=$(wc -l < {input.stats}/stats.tsv) '{{if(NR==(total-2)/2){{print $6}}; \
            if(NR==(total-1)/2){{print $6;exit}}}}') > {log} 2>&1

        echo -e "TargetDepth: $TargetDepth" >> {log}
        qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --i-phylogeny {input.tree} \
            --p-max-depth $TargetDepth \
            --p-min-depth 1 \
            --p-steps {params.p_steps} \
            --p-iterations {params.p_iters} \
            {params.metadata} \
            --o-visualization {output} >> {log} 2>&1
        """


use rule demux_summary_export as alpha_rarefaction_export with:
    input:
        rules.alpha_rarefaction.output
    output:
        directory('diversity/alpha-rarefaction')
    benchmark:
        '.log/diversity/alpha_rarefaction_export.bm'
    log:
        '.log/diversity/alpha_rarefaction_export.log'
