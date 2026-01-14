
process run_parameters {
  // This process generates a parameters file that includes the operational
  // parameters used in the workflow.

    tag '12'
    publishDir params.outDir + '/parameters', mode: 'copy'

    output:
        path("parameters.txt")

    script:
        """
        touch parameters.tsv
        echo "Parameter\tValue" >> parameters.txt
        echo "Mapping quality:\t>=${params.mapqual}" >> parameters.txt
        echo "Mutation frequency:\t>=${params.varthres}" >> parameters.txt
        echo "Variant quality:\t>${params.varqual}" >> parameters.txt
        echo "Minimum sequence depth:\t${params.depth}" >> parameters.txt
        echo "BWA seed length:\t${params.seed_length}" >> parameters.txt
        """
}