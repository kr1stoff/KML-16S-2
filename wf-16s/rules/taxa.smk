rule taxa:
    input:
        rules.dada2_rename_rep_seqs.output.re_rep_seqs
    output:
        taxonomy='taxa/taxonomy.qza',
        taxonomy_dir=directory('taxa/taxonomy')
    benchmark:
        '.log/taxa/taxa.bm'
    log:
        '.log/taxa/taxa.log'
    conda:
        config['conda']['qiime2']
    params:
        # 内存不够, 需要修改 --p-reads-per-batch 参数, 默认 20000
        # 当前参数 10 min
        '--p-confidence 0.7 --p-reads-per-batch 15000'
    threads:
        config['threads']['high']
    resources:
        tmpdir='.temp/qiime2'
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {config[database][silva]} \
            --i-reads {input} \
            --p-n-jobs {threads} \
            {params} \
            --o-classification {output.taxonomy} > {log} 2>&1 && \
        qiime tools export \
            --input-path {output.taxonomy} \
            --output-path {output.taxonomy_dir} >> {log} 2>&1
        """


rule taxa_visual:
    input:
        rules.taxa.output.taxonomy
    output:
        'taxa/taxonomy.qzv'
    benchmark:
        '.log/taxa/taxa_visual.bm'
    log:
        '.log/taxa/taxa_visual.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime metadata tabulate \
            --m-input-file {input} \
            --o-visualization {output} > {log} 2>&1
        """


use rule demux_summary_export as taxa_visual_export with:
    input:
        rules.taxa_visual.output
    output:
        directory('taxa/taxonomy_v')
    benchmark:
        '.log/taxa/taxa_visual_export.bm'
    log:
        '.log/taxa/taxa_visual_export.log'


rule taxa_barplot:
    input:
        table=rules.dada2_rename_table.output.re_table,
        taxonomy=rules.taxa.output.taxonomy
    output:
        'taxa/taxonomy_barplot.qzv'
    benchmark:
        '.log/taxa/taxa_barplot.bm'
    log:
        '.log/taxa/taxa_barplot.log'
    conda:
        config['conda']['qiime2']
    params:
        '--m-metadata-file ' + config['metadata']
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            {params} \
            --o-visualization {output} > {log} 2>&1
        """


use rule demux_summary_export as taxa_barplot_export with:
    input:
        rules.taxa_barplot.output
    output:
        directory('taxa/taxonomy_barplot_v')
    benchmark:
        '.log/taxa/taxa_barplot_export.bm'
    log:
        '.log/taxa/taxa_barplot_export.log'


from enum import Enum


class RanksLevel(Enum):
    Kingdom = 1
    Phylum = 2
    Class = 3
    Order = 4
    Family = 5
    Genus = 6
    Species = 7


rule taxa_collapse:
    input:
        table=rules.diversity_core.output.rarefied_table,
        taxonomy=rules.taxa.output.taxonomy
    output:
        table='taxa/collapse/{levels}.qza',
        table_dir=directory('taxa/collapse/{levels}'),
        table_tsv='taxa/collapse/{levels}/feature-table.tsv'
    benchmark:
        '.log/taxa/taxa_collapse_{levels}.bm'
    log:
        '.log/taxa/taxa_collapse_{levels}.log'
    conda:
        config['conda']['qiime2']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    params:
        levels=lambda wildcards: RanksLevel[wildcards.levels].value
    shell:
        """
        qiime taxa collapse \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-level {params.levels} \
            --o-collapsed-table {output.table} > {log} 2>&1 && \
        qiime tools export \
            --input-path {output.table} \
            --output-path {output.table_dir} >> {log} 2>&1 && \
        biom convert \
            -i {output.table_dir}/feature-table.biom \
            -o {output.table_tsv} --to-tsv >> {log} 2>&1
        """


rule taxa_collapse_freq:
    input:
        rules.taxa_collapse.output.table
    output:
        table='taxa/collapse_freq/{levels}.qza',
        table_dir=directory('taxa/collapse_freq/{levels}'),
        table_tsv='taxa/collapse_freq/{levels}/feature-table.tsv'
    benchmark:
        '.log/taxa/taxa_collapse_freq_{levels}.bm'
    log:
        '.log/taxa/taxa_collapse_freq_{levels}.log'
    conda:
        config['conda']['qiime2']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    shell:
        """
        qiime feature-table relative-frequency \
            --i-table {input} \
            --o-relative-frequency-table {output.table} > {log} 2>&1 && \
        qiime tools export \
            --input-path {output.table} \
            --output-path {output.table_dir} >> {log} 2>&1 && \
        biom convert \
            -i {output.table_dir}/feature-table.biom \
            --header-key "taxonomy" \
            -o {output.table_tsv} --to-tsv >> {log} 2>&1
        """


rule taxa_barplot_R:
    input:
        rules.taxa_collapse.output.table_tsv
    output:
        img_dir=directory('taxa/barplot/{levels}'),
        result='taxa/barplot/{levels}/{levels}.result.txt',
        barstat='taxa/barplot/{levels}/{levels}.barstat.txt',
        img_pdf='taxa/barplot/{levels}/{levels}.barplot.pdf',
        img_png='taxa/barplot/{levels}/{levels}.barplot.png'
    benchmark:
        '.log/taxa/taxa_barplot_R_{levels}.bm'
    log:
        '.log/taxa/taxa_barplot_R_{levels}.log'
    conda:
        config['conda']['microplot']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    shell:
        """
        mkdir -p {output.img_dir}
        Rscript {config[my_scripts]}/taxa_barplot.R \
            {input} \
            {output.img_dir}/{wildcards.levels} > {log} 2>&1
        """


rule krona_plot:
    input:
        rules.taxa_collapse.output.table_tsv
    output:
        krona_input='taxa/krona/{sample}.krona.{levels}.txt',
        krona_html='taxa/krona/{sample}.krona.{levels}.html',
        krona_pdf='taxa/krona/{sample}.krona.{levels}.pdf',
        krona_png='taxa/krona/{sample}.krona.{levels}.png'
    resources:
        time='1h',# 设置运行的最大时间限制
        mem_mb='1000'  # 设置最大内存使用限制(Mb)
    benchmark:
        '.log/taxa/{sample}_krona_plot_{levels}.bm'
    log:
        '.log/taxa/{sample}_krona_plot_{levels}.log'
    conda:
        config['conda']['microplot']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    shell:
        """
        python {config[my_scripts]}/table_to_krona_input.py {input} {wildcards.sample} > {output.krona_input} 2> {log}
        ktImportText {output.krona_input} -o {output.krona_html} >> {log} 2>&1
        export OPENSSL_CONF=/dev/null
        phantomjs {config[my_scripts]}/snapshot_krona_html.js {output.krona_html} {output.krona_pdf} >> {log} 2>&1
        phantomjs {config[my_scripts]}/snapshot_krona_html.js {output.krona_html} {output.krona_png} >> {log} 2>&1
        """


rule heatmap_hclust2_plot:
    input:
        rules.taxa_collapse_freq.output.table_tsv
    output:
        hclust2_labelsize='taxa/heatmap_hclust2/{levels}.hclust2.labelsize',
        hclust2_input='taxa/heatmap_hclust2/{levels}.hclust2.txt',
        hclust2_img='taxa/heatmap_hclust2/{levels}.hclust2.png'
    benchmark:
        '.log/taxa/heatmap_hclust2_plot_{levels}.bm'
    log:
        '.log/taxa/heatmap_hclust2_plot_{levels}.log'
    conda:
        config['conda']['graphlan']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    params:
        f_dist_f='braycurtis',#特征值绘制树/聚类方法
        s_dist_f='braycurtis',#样本绘制树/聚类方法
        #cell_aspect_ratio = '0.5', #每个小的矩形单云色块的 宽/高 比值
        log_scale='--log_scale',#是否取对数值绘图, 取对数可以把 过大/过小 的值很好地显示出来
        #flabel_size = 2,      #特征物种名称标签大小
        #slabel_size = 2,      #样本标签大小
        max_flabel_len=100,#物种名称标签显示的最长字符数
        max_slabel_len=30,#样本标签显示的最长字符数
        minv=0.001,#数据最低值显示为NaN颜色
        dpi=300,#dpi像素
        cmode='YlGn',#色卡模板
        ftop='100'  #选取特征值的前ftop个绘图
    shell:
        """
        python {config[my_scripts]}/table_to_hclust2_input.py \
            {input} {params.ftop} > {output.hclust2_input} 2> {output.hclust2_labelsize}
        flabelsize=$(cat {output.hclust2_labelsize} | cut -f 2)
        slabelsize=$(cat {output.hclust2_labelsize} | cut -f 1)
        cellratio=$(cat {output.hclust2_labelsize} | cut -f 3)
        hclust2.py -i {output.hclust2_input} \
            -o {output.hclust2_img} \
            --f_dist_f {params.f_dist_f} \
            --s_dist_f {params.s_dist_f} \
            --cell_aspect_ratio $cellratio \
            --flabel_size $flabelsize \
            --slabel_size $slabelsize \
            --max_flabel_len {params.max_flabel_len} \
            --max_slabel_len {params.max_slabel_len} \
            --minv {params.minv} \
            --dpi {params.dpi} \
            -c {params.cmode} \
            --ftop {params.ftop} \
            {params.log_scale} > {log} 2>&1
        """


rule graphlan_plot:
    input:
        rules.taxa_collapse_freq.output.table_tsv
    output:
        graphlan_input='taxa/graphlan/{levels}.graphlan.txt',
        graphlan_tree='taxa/graphlan/{levels}.graphlan.tree.txt',
        graphlan_annot='taxa/graphlan/{levels}.graphlan.annot.txt',
        graphlan_xml='taxa/graphlan/{levels}.graphlan.xml',
        graphlan_img='taxa/graphlan/{levels}.graphlan.png'
    benchmark:
        '.log/taxa/hgraphlan_plot_{levels}.bm'
    log:
        '.log/taxa/hgraphlan_plot_{levels}.log'
    conda:
        config['conda']['graphlan']
    wildcard_constraints:
        levels='(Kingdom)|(Phylum)|(Class)|(Order)|(Family)|(Genus)|(Species)'
    shell:
        """
        python {config[my_scripts]}/table_to_graphlan_input.py \
            {input} > {output.graphlan_input} 2> {log}
        export2graphlan.py -i {output.graphlan_input} \
            --tree {output.graphlan_tree} \
            --annotation {output.graphlan_annot} \
            --most_abundant 100 \
            --abundance_threshold 1 \
            --least_biomarkers 10 \
            --annotations 5,6 \
            --external_annotations 7 \
            --min_clade_size 1 > {log} 2>&1
        graphlan_annotate.py \
            --annot {output.graphlan_annot} \
            {output.graphlan_tree} {output.graphlan_xml} > {log} 2>&1
        graphlan.py {output.graphlan_xml} {output.graphlan_img} --dpi 300 --external_legends > {log} 2>&1
        """
