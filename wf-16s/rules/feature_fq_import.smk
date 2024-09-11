from pathlib import Path


rule generate_manifest:
    input:
        fq1=expand('qc/fastp/{sample}.1.fastq.gz',sample=config['samples']),
        fq2=expand('qc/fastp/{sample}.2.fastq.gz',sample=config['samples'])
    output:
        'feature/clean_fq.manifest'
    run:
        with open(output[0],'w') as f:
            f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")

            for fq1, fq2 in zip(input.fq1,input.fq2):
                sample = Path(fq1).stem.split('.')[0]
                f.write(f'{sample}\t{Path(fq1).resolve()}\t{Path(fq2).resolve()}\n')


rule qiime2_fq_import:
    input:
        rules.generate_manifest.output
    output:
        demux_qza='feature/demux.qza',
        demux_qzv='feature/demux.qzv'
    benchmark:
        '.log/feature/qiime2_fq_import.bm'
    log:
        '.log/feature/qiime2_fq_import.log'
    params:
        "--type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2"
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime tools import \
            {params} \
            --input-path {input} \
            --output-path {output.demux_qza} \
            >> {log} 2>&1
        qiime demux summarize \
            --i-data {output.demux_qza} \
            --o-visualization {output.demux_qzv} >> {log} 2>&1
        """


rule demux_summary_export:
    input:
        rules.qiime2_fq_import.output.demux_qzv
    output:
        directory('feature/demux_v')
    benchmark:
        '.log/feature/demux_summary_export.bm'
    log:
        '.log/feature/demux_summary_export.log'
    conda:
        config['conda']['qiime2']
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output} > {log} 2>&1
        """
