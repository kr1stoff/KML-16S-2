rule picrust2_pipeline:
    input:
        seqs_fa=rules.dada2_rename_rep_seqs.output.re_rep_seqs_fa,
        otu_tsv=rules.diversity_core_export.output
    output:
        fdir=directory('function'),
        pdir=directory('function/picrust2_out'),
        pathway_tsv='function/picrust2_out/path_abun_unstrat.tsv'
    benchmark:
        '.log/diff/picrust2_pipeline.bm'
    log:
        '.log/diff/picrust2_pipeline.log'
    threads:
        config['threads']['low']
    conda:
        config['conda']['qiime2']
    shell:
        """
        grep -v "^# Constructed from biom file" {input.otu_tsv} > {output.fdir}/picrust2.otu.input 2> {log}
        if [ -d "{output.pdir}" ];then
            rm -rf {output.pdir}
        fi
        picrust2_pipeline.py -s {input.seqs_fa} -i {output.fdir}/picrust2.otu.input \
            -o {output.pdir} -p {threads} >> {log} 2>&1
        gunzip -dc {output.pdir}/pathways_out/path_abun_unstrat.tsv.gz > {output.pathway_tsv} 2>> {log}
        """


rule picrust2_plot:
    input:
        pathway_tsv=rules.picrust2_pipeline.output.pathway_tsv,
        metadata=rules.select_sampling_depth.output.re_metadata
    output:
        directory('function/picrust2_plot')
    benchmark:
        '.log/diff/picrust2_plot.bm'
    log:
        '.log/diff/picrust2_plot.log'
    conda:
        config['conda']['microplot']
    shell:
        """
        mkdir -p {output}
        Rscript {config[my_scripts]}/ggpicrust2.R {input.metadata} {input.pathway_tsv} group {output}/group > {log} 2>&1
        """
