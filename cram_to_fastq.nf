#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CRAM2BAM {
    input:
    tuple val(sample_id), file(cram)
    val(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_ch

    script:
    """
    samtools view -b -T $reference -o ${sample_id}".bam"  $cram    
    """
}

process BAM2FASTQ {
    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq"), path("${sample_id}_2.fastq"), emit: fastq_ch

    script:
    """
    bamtofastq collate=1 inputformat=bam exclude=SECONDARY,SUPPLEMENTARY \
            F=${sample_id}_1.fastq \
            O=${sample_id}_1_orphan.fastq \
            F2=${sample_id}_2.fastq \
            O2=${sample_id}_2_orphan.fastq \
            S=${sample_id}_single_end.fastq \
            < ${bam}
    """
}

process GZIP_FASTQ {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    input:
    tuple val(sample_id), file(read_1), file(read_2)

    output:
    path("*.fastq.gz")

    script:
    """
    gzip -f *.fastq
    """
}

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    cram_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.cram)) }

    CRAM2BAM(cram_path_ch, params.reference)

    BAM2FASTQ(CRAM2BAM.out.bam_ch)

    GZIP_FASTQ(BAM2FASTQ.out.fastq_ch)
}    
