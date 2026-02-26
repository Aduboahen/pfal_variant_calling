process mark_duplicates {
// mark duplicate sequences with picard

    tag '3'
    publishDir params.outDir + "/mapping", mode: 'copy'

    input:
        path(bam)

    output:
        path("${params.sampleid}.markdup.bam", emit: bam)
        path("${params.sampleid}.markdup.bam.bai", emit: bam_index)
        path("${params.sampleid}.markdup.stats")

    script:
    """
    picard MarkDuplicates -I ${bam} -O ${params.sampleid}.markdup.bam -M ${params.sampleid}.metrics.txt
    samtools index ${params.sampleid}.markdup.bam > ${params.sampleid}.markdup.bam.bai
    samtools stats ${params.sampleid}.markdup.bam > ${params.sampleid}.markdup.stats
    """
}