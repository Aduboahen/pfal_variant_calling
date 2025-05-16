process mapping {
    // map reads to reference and output sorted bam files

    tag '2'
    publishDir params.outDir + '/mapping', mode: 'copy'

    input:
      path(trimmed1)
      path(trimmed2)
      path(merged)

    output:
      path("${params.sampleid}.bam"), emit: bam
      path("${params.sampleid}.bam.bai"), emit: bam_index
      path("${params.sampleid}.stats")

    script:

      """
        bwa mem -M -k 10 -t ${params.threads} ${params.reference} ${merged} | samtools sort -T temp -O bam -o ${params.sampleid}.long.bam

        bwa mem -M -k 10 -t ${params.threads} ${params.reference} ${trimmed1} ${trimmed2} | samtools sort -T temp -O bam -o ${params.sampleid}.short.bam

        samtools merge ${params.sampleid}.merged.bam ${params.sampleid}.long.bam ${params.sampleid}.short.bam

        samtools sort ${params.sampleid}.merged.bam -@ ${params.threads} -o ${params.sampleid}.bam

        samtools index ${params.sampleid}.bam > ${params.sampleid}.bam.bai

        samtools stats ${params.sampleid}.bam > ${params.sampleid}.stats
      """
}