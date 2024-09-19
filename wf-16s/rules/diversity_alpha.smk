rule diversity_alpha:
    input:
        rules.diversity_core.output.rarefied_table
    output:
        metric=expand('diversity/alpha/pmetric.{p_metric}.qza',
            p_metric=['ace', 'chao1', 'simpson', 'shannon', 'goods_coverage', 'observed_features']),
        metric_dir=directory(expand('diversity/alpha/pmetric.{p_metric}',
            p_metric=['ace', 'chao1', 'simpson', 'shannon', 'goods_coverage', 'observed_features']))
    benchmark:
        '.log/diversity/diversity_alpha.bm'
    log:
        '.log/diversity/diversity_alpha.log'
    conda:
        config['conda']['qiime2']
    params:
        p_metric=['ace', 'chao1', 'simpson', 'shannon', 'goods_coverage', 'observed_features'],
        pm_outdir='diversity/alpha'
    shell:
        """
        echo -n > {log}
        for pmetric in {params.p_metric}
        do
        qiime diversity alpha \
            --i-table {input} \
            --o-alpha-diversity {params.pm_outdir}/pmetric.$pmetric.qza \
            --p-metric $pmetric >> {log} 2>&1
        qiime tools export \
            --input-path {params.pm_outdir}/pmetric.$pmetric.qza \
            --output-path {params.pm_outdir}/pmetric.$pmetric >> {log} 2>&1
        done
        """


rule diversity_alpha_stats:
    input:
        rules.diversity_alpha.output.metric_dir
    output:
        'diversity/alpha_diversity_stats.tsv'
    benchmark:
        '.log/diversity/diversity_alpha_stats.bm'
    log:
        '.log/diversity/diversity_alpha_stats.log'
    conda:
        config['conda']['qiime2']
    params:
        p_metrics=['ace', 'chao1', 'simpson', 'shannon', 'goods_coverage', 'observed_features'],
        pm_outdir='diversity/alpha'
    shell:
        """
        tsvs=""
        for pmetric in {params.p_metrics}
        do
            tsvs=$tsvs\" \"{params.pm_outdir}/pmetric.$pmetric/alpha-diversity.tsv
        done
        python {config[my_scripts]}/combine_alpha_diversity_stat.py $tsvs > {output} 2> {log}
        """
