rule diversity_beta:
    input:
        rules.diversity_core.output.rarefied_table
    output:
        metric=expand('diversity/beta/pmetric.{p_metric}.qza',p_metric=['braycurtis', 'jaccard', 'euclidean']),
        metric_dir=directory(expand('diversity/beta/pmetric.{p_metric}',
            p_metric=['braycurtis', 'jaccard', 'euclidean']))
    benchmark:
        '.log/diversity/diversity_beta.bm'
    log:
        '.log/diversity/diversity_beta.log'
    threads:
        config['threads']['low']
    conda:
        config['conda']['qiime2']
    params:
        tools='beta',
        p_metrics=['braycurtis', 'jaccard', 'euclidean'],
        pm_outdir='diversity/beta',
        n_jobs='--p-n-jobs 4'
    shell:
        """
        echo -n > {log}
        for pmetric in {params.p_metrics}
        do
        qiime diversity beta \
            --i-table {input} \
            --o-distance-matrix {params.pm_outdir}/pmetric.$pmetric.qza \
            {params.n_jobs} \
            --p-metric $pmetric >> {log} 2>&1
        qiime tools export \
            --input-path {params.pm_outdir}/pmetric.$pmetric.qza \
            --output-path {params.pm_outdir}/pmetric.$pmetric >> {log} 2>&1
        done
        """
