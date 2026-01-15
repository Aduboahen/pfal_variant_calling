process filter {
  // based on params.depth, quality, mapping quality and variant threshold

  tag '6'
  publishDir params.outDir + "/variants", mode: 'copy'

  input:
  path vcf

  output:
  path ("${params.sampleid}.vcf"), emit: snps

  script:
  """
      bcftools filter -i 'QUAL>=${params.varqual} \
      && FORMAT/DP>${params.depth} \
      && MQ>=${params.mapqual} \
      && DP4[2] > 0 \
      && DP4[3] > 0 \
      && DP4[2]/(DP4[2]+DP4[0])>=${params.varthres}\
      && DP4[3]/(DP4[3]+DP4[1])>=${params.varthres}'\
      -g10 \
      -o ${params.sampleid}.filtered.vcf ${vcf}


      # zcat ${params.sampleid}.filtered.vcf | sed -e 's/Pf3D7_01_v3/1/g;\
      # s/Pf3D7_02_v3/2/g; \
      # s/Pf3D7_03_v3/3/g; \
      # s/Pf3D7_04_v3/4/g; \
      # s/Pf3D7_05_v3/5/g; \
      # s/Pf3D7_06_v3/6/g; \
      # s/Pf3D7_07_v3/7/g; \
      # s/Pf3D7_08_v3/8/g; \
      # s/Pf3D7_09_v3/9/g; \
      # s/Pf3D7_10_v3/10/g; \
      # s/Pf3D7_11_v3/11/g; \
      # s/Pf3D7_12_v3/12/g; \
      # s/Pf3D7_13_v3/13/g; \
      # s/Pf3D7_14_v3/14/g; \
			# /Pf3D7_API_v3/d; \
      # /Pf3D7_MIT_v3/d' | bgzip > ${params.sampleid}.vcf

      bcftools index -f ${params.sampleid}.vcf -o ${params.sampleid}.vcf.csi
    """
}

process filter_snps {
  // filter snps
  // based on params.depth, quality, mapping quality and variant threshold

  tag '6'
  publishDir params.outDir + "/variants", mode: 'copy'

  input:
  path vcf

  output:
  path ("${params.sampleid}.bcftools.vcf.gz")

  script:
  """
      bcftools filter -i 'type="snp" \
      && QUAL>=${params.varqual} \
      && FORMAT/DP>${params.depth} \
      && DP4[2] > 0 \
      && DP4[3] > 0 \
      && DP4[2]/(DP4[2]+DP4[0])>=${params.varthres} \
      && DP4[3]/(DP4[3]+DP4[1])>=${params.varthres}' \
      -Oz \
      -o ${params.sampleid}.bcftools.vcf.gz ${vcf}

      bcftools index -f ${params.sampleid}.bcftools.vcf.gz
    """
}

process filter_indels {
  // filter indels 

  tag '6b'
  publishDir params.outDir + "/variants", mode: 'copy'

  input:
  path vcf

  output:
  path ("${params.sampleid}.indels.vcf"), emit: indels
  path ("${params.sampleid}.indels.vcf.csi"), emit: indels_index

  script:
  """
      bcftools filter -i 'type="indel"\
      && QUAL>=${params.varqual} \
      && FORMAT/DP>${params.depth} \
      && DP4[2] > 0 \
      && DP4[3] > 0 \
      && DP4[2]/(DP4[2]+DP4[0])>=${params.varthres}\
      && DP4[3]/(DP4[3]+DP4[1])>=${params.varthres}'\
      -g10 \
      -O z \
      -o ${params.sampleid}.indels.filtered.vcf  ${vcf}

      bcftools index -f ${params.sampleid}.indels.vcf -o ${params.sampleid}.indels.vcf.csi
    """
}

process merge_vcf {
  // merge snp and indel vcf files

  tag '6c'
  publishDir params.outDir + "/variants", mode: 'copy'

  input:
  path snps
  path indels
  path snps_index
  path indels_index

  output:
  path "${params.sampleid}.merged.vcf"

  script:
  """
      bcftools merge --threads ${params.threads} -O z -o ${params.sampleid}.merged.vcf ${snps} ${indels}
    """
}
process filter_lofreq {
  // based on params.depth, quality, mapping quality and variant threshold

  tag '6'
  publishDir params.outDir + "/variants", mode: 'copy'

  input:
  path vcf

  output:
  path "${params.sampleid}.lofreq.filtered.vcf.gz"

  script:
  """
      lofreq filter -i ${vcf} --only-snvs --print-all \
      --cov-min ${params.depth} \
      --af-min ${params.varthres} \
      --snvqual-thresh ${params.varqual} \
      --sb-thresh ${params.strandthres} \
      -o ${params.sampleid}.lofreq.filtered.vcf.gz

      bcftools index -f ${params.sampleid}.lofreq.filtered.vcf.gz
    """
}
