rule phylogeny:
    input:
        rules.dada2_rename_rep_seqs.output.re_rep_seqs
    output:
        aligned_rep_seqs='feature/phylogeny/aligned-rep-seqs.qza',
        masked_aligned_rep_seqs='feature/phylogeny/masked-aligned-rep-seqs.qza',
        unrooted_tree='feature/phylogeny/unrooted-tree.qza',
        rooted_tree='feature/phylogeny/rooted-tree.qza'
    benchmark:
        '.log/phylogeny/phylogeny.bm'
    log:
        '.log/phylogeny/phylogeny.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input} \
            --o-alignment {output.aligned_rep_seqs}\
            --o-masked-alignment {output.masked_aligned_rep_seqs} \
            --o-tree {output.unrooted_tree} \
            --o-rooted-tree {output.rooted_tree} > {log} 2>&1
        """


use rule demux_summary_export as phylogeny_aligned_repseqs_export with:
    input:
        rules.phylogeny.output.aligned_rep_seqs
    output:
        directory('feature/phylogeny/aligned-rep-seqs')
    benchmark:
        '.log/phylogeny/phylogeny_aligned_repseqs_export.bm'
    log:
        '.log/phylogeny/phylogeny_aligned_repseqs_export.log'


use rule demux_summary_export as phylogeny_masked_aligned_repseqs_export with:
    input:
        rules.phylogeny.output.masked_aligned_rep_seqs
    output:
        directory('feature/phylogeny/masked-aligned-rep-seqs')
    benchmark:
        '.log/phylogeny/phylogeny_masked_aligned_repseqs_export.bm'
    log:
        '.log/phylogeny/phylogeny_masked_aligned_repseqs_export.log'


use rule demux_summary_export as phylogeny_unrooted_tree_export with:
    input:
        rules.phylogeny.output.unrooted_tree
    output:
        directory('feature/phylogeny/unrooted-tree')
    benchmark:
        '.log/phylogeny/phylogeny_unrooted_tree_export.bm'
    log:
        '.log/phylogeny/phylogeny_unrooted_tree_export.log'


use rule demux_summary_export as phylogeny_rooted_tree_export with:
    input:
        rules.phylogeny.output.rooted_tree
    output:
        directory('feature/phylogeny/rooted-tree')
    benchmark:
        '.log/phylogeny/phylogeny_rooted_tree_export.bm'
    log:
        '.log/phylogeny/phylogeny_rooted_tree_export.log'


rule phylogeny_treeplot:
    input:
        rules.phylogeny_rooted_tree_export.output
    output:
        'feature/phylogeny/rooted-tree.{layout}.{suffix}'
    benchmark:
        '.log/phylogeny/rooted-tree.{layout}.{suffix}.bm'
    log:
        '.log/phylogeny/rooted-tree.{layout}.{suffix}.log'
    wildcard_constraints:
        layout='(rectangular)|(circular)',
        suffix='(pdf)|(png)'
    conda:
        config['conda']['microplot']
    params:
        '--branchlength none'  #none, branch.length
    shell:
        """
        Rscript {config[my_scripts]}/phylogeny_tree_plot.R \
            -i {input}/tree.nwk \
            -o {output} \
            --layout {wildcards.layout} \
            {params} > {log} 2>&1
        """
