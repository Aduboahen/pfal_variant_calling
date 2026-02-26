process genome_stats {
    // This process will generate a file with the following stats:
    // 1. Number of sequences
    // 2. Total length
    // 3. GC content
    // Genome stats with faCount
    // 4. Number of N

    tag '11'
    publishDir params.outDir + '/consensus', mode: 'copy'

    input:
        path(fasta)

    output:
        path("${params.sampleid}.stats.txt")

    script:
        """
        faCount ${fasta} > ${params.sampleid}.stats.txt
        ${params.parse_stats} --stats ${params.sampleid}.stats.txt
        """
}
