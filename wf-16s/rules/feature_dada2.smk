rule qiime2_dada2:
    input:
        rules.qiime2_fq_import.output.demux_qza
    output:
        stats='dada2/stats.qza',
        table='dada2/raw-table.qza',
        rep_seqs='dada2/raw-rep-seqs.qza'
    benchmark:
        '.log/dada2/qiime2_dada2.bm'
    log:
        '.log/dada2/qiime2_dada2.log'
    conda:
        config['conda']['qiime2']
    threads:
        config['threads']['high']
    params:
        '--p-trim-left 0 --p-trunc-len 0'
    shell:
        """
        qiime dada2 denoise-single \
            {params} \
            --i-demultiplexed-seqs {input} \
            --p-n-threads {threads} \
            --o-representative-sequences {output.rep_seqs} \
            --o-denoising-stats {output.stats} \
            --o-table {output.table} > {log} 2>&1
        """


rule dada2_rename_table:
    input:
        rules.qiime2_dada2.output.table
    output:
        table_dir=directory('dada2/raw-table'),
        tsv='dada2/raw-table/feature-table.tsv',
        re_tsv='dada2/table/feature-table.tsv',
        re_biom='dada2/raw-table/rename.feature-table.biom',
        re_table='dada2/table.qza'
    benchmark:
        '.log/dada2/dada2_rename_table.bm'
    log:
        '.log/dada2/dada2_rename_table.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output.table_dir} > {log} 2>&1 
        biom convert -i {output.table_dir}/feature-table.biom \
            -o {output.tsv} --to-tsv >> {log} 2>&1 
        cat {output.tsv} | awk '{{
            if($0~/#/)
            {{
                print $0
            }}
            else{{
                seqid+=1;
                $1="ASV_"seqid;
                gsub(" ", "\\t", $0);
                print $0
            }}
            }}' > {output.re_tsv} 2>> {log}
        biom convert -i {output.re_tsv} -o {output.re_biom} --to-hdf5 >> {log} 2>&1
        qiime tools import \
            --input-path {output.re_biom} \
            --type "FeatureTable[Frequency]" \
            --output-path {output.re_table} >> {log} 2>&1
        """


rule dada2_table_freq:
    input:
        rules.dada2_rename_table.output.re_table
    output:
        table_freq="dada2/table_freq.qza",
        table_freq_dir=directory("dada2/table_freq"),
        table_freq_tsv="dada2/table_freq/feature-table.tsv"
    benchmark:
        '.log/dada2/dada2_table_freq.bm'
    log:
        '.log/dada2/dada2_table_freq.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime feature-table relative-frequency \
            --i-table {input} \
            --o-relative-frequency-table {output.table_freq} > {log} 2>&1 && \
        qiime tools export \
            --input-path {output.table_freq} \
            --output-path {output.table_freq_dir} >> {log} 2>&1 && \
        biom convert \
            -i {output.table_freq_dir}/feature-table.biom \
            -o {output.table_freq_tsv} \
            --to-tsv >> {log} 2>&1         
        """


rule dada2_rename_rep_seqs:
    input:
        rules.qiime2_dada2.output.rep_seqs
    output:
        rep_seqs_dir=directory('dada2/raw-rep-seqs'),
        re_rep_seqs_fa='dada2/rep-seqs/dna-sequences.fasta',
        re_rep_seqs='dada2/rep-seqs.qza'
    benchmark:
        '.log/dada2/dada2_rename_rep_seqs.bm'
    log:
        '.log/dada2/dada2_rename_rep_seqs.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output.rep_seqs_dir} > {log} 2>&1
        cat {output.rep_seqs_dir}/dna-sequences.fasta | awk '{{
            if($0~/^>/)
            {{
                seqid+=1;
                print ">ASV_"seqid;
            }}
            else{{
                print $0
            }}
            }}' > {output.re_rep_seqs_fa} 2>> {log}
        qiime tools import \
            --input-path {output.re_rep_seqs_fa} \
            --type "FeatureData[Sequence]" \
            --output-path {output.re_rep_seqs} >> {log} 2>&1
        """


use rule demux_summary_export as stats_export with:
    input:
        rules.qiime2_dada2.output.stats
    output:
        directory('dada2/stats')
    benchmark:
        '.log/dada2/stats_export.bm'
    log:
        '.log/dada2/stats_export.log'


rule stats_visual:
    input:
        rules.qiime2_dada2.output.stats
    output:
        'dada2/stats.qzv'
    benchmark:
        '.log/dada2/stats_visual.bm'
    log:
        '.log/dada2/stats_visual.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime metadata tabulate \
            --m-input-file {input} \
            --o-visualization {output} > {log} 2>&1
        """


use rule demux_summary_export as stats_visual_export with:
    input:
        rules.stats_visual.output
    output:
        directory('dada2/stats_v')
    benchmark:
        '.log/dada2/stats_visual_export.bm'
    log:
        '.log/dada2/stats_visual_export.log'


rule rep_seqs_visual:
    input:
        rules.dada2_rename_rep_seqs.output.re_rep_seqs
    output:
        'dada2/rep-seqs.qzv'
    benchmark:
        '.log/dada2/rep_seqs_visual.bm'
    log:
        '.log/dada2/rep_seqs_visual.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime feature-table tabulate-seqs \
            --i-data {input} \
            --o-visualization {output} > {log} 2>&1
        """


use rule demux_summary_export as rep_seqs_visual_export with:
    input:
        rules.rep_seqs_visual.output
    output:
        directory('dada2/rep-seqs_v')
    benchmark:
        '.log/dada2/rep_seqs_visual_export.bm'
    log:
        '.log/dada2/rep_seqs_visual_export.log'


rule table_visual:
    input:
        rules.dada2_rename_table.output.re_table
    output:
        'dada2/table.qzv'
    benchmark:
        '.log/dada2/table_visual.bm'
    log:
        '.log/dada2/table_visual.log'
    conda:
        config['conda']['qiime2']
    params:
        '--m-sample-metadata-file ' + config['metadata']
    shell:
        """
        qiime feature-table summarize \
            --i-table {input} \
            --o-visualization {output} \
            {params} > {log} 2>&1
        """


use rule demux_summary_export as table_visual_export with:
    input:
        rules.table_visual.output
    output:
        directory('dada2/table_v')
    benchmark:
        '.log/dada2/table_visual_export.bm'
    log:
        '.log/dada2/table_visual_export.log'
