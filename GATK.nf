#!/usr/bin/env nextflow

nextflow.enable.dsl=1

// parameters to pass to command

outDir           = "$params.outDir"
params.reads     = "$projectDir/test_fq/test_{1,2}.fastq.gz"
reference        = "$projectDir/reference/Pf3D7.fasta"
known_sites      = "$projectDir/reference/3d7_hb3.gatk.final.vcf"
sampleid         = "$params.sampleid"
params.threads   = 6 
threads          = "$params.threads"


// tools

vcf2table        = "$projectDir/scripts/vcf2table_gatk.py"
parse_stats      = "$projectDir/scripts/parse_stats.py"

// filters

params.QD         = 2.0
QD                = "$params.QD" 
params.FS         = 60.0
FS                = "$params.FS"
params.MQRankSum  = -12.5 
MQRankSum         = "$params.MQRankSum" 
params.MMQ        =  15.0
MMQ               = "$params.MMQ" 
params.DP         = 5
params.MQ         =  40.0
MQ               = "$params.MQ" 
params.DP         = 5
DP                = "$params.DP"
SOR               = 3.0
ReadPosRankSum    = -8.0

// channel to get reads as tuples

Channel
    .fromFilePairs(params.reads)
    .ifEmpty{"no such files"}
    .set{reads_ch}

log.info"""
P. falciparum variant calling

===============Sources===================
codeBase            : "$projectDir"
Sample              : "$params.sampleid"
Reads               : "$params.reads"
Output Directory    : "$params.outDir"
Threads             : "$params.threads"


==============Variant filters============
Coverage                : "$params.DP"
Quality By Depth        : "$params.QD" 
Fisher Strand Bias      : "$params.FS"
Median Mapping Quality  : "$params.MMQ" 
RMS Mapping Quality     : "$params.MQ" 


==============Reference Sequence=========
Reference : "Pf3D7"


==============Author======================
James Osei-Mensa
oseimensa@kccr.de
"""

// Clean reads (adapter and read length filter)

process 'clean reads' {
    tag '1A'
    publishDir outDir + '/QC', mode: 'copy'
    input:
      set pairID, file(reads) from reads_ch
    output:
      set file("${reads[0].baseName}_fastp.fastq.gz"), file("${reads[1].baseName}_fastp.fastq.gz") into fastp_2A
      file "${sampleid}_fastpmerged.fastq.gz" into fastpmerged_2A
      file "${sampleid}.fastp.json"
    script:
    """
      fastp -i ${reads[0]} -I ${reads[1]} \
      -o ${reads[0].baseName}_fastp.fastq.gz -O ${reads[1].baseName}_fastp.fastq.gz \
      --merge --merged_out ${sampleid}_fastpmerged.fastq.gz \
      --overlap_diff_limit 0 \
      --trim_poly_x --trim_poly_g --length_required 100 --thread ${threads} --detect_adapter_for_pe\
      --json ${sampleid}.fastp.json --html ${sampleid}.fastp.html
    """
}

// map reads to reference and output sorted bam files

process mapping {
    tag '2A'
    publishDir outDir + "/mapping", mode: 'copy'
    input:
      set file(trimmed1), file(trimmed2) from fastp_2A
      file fastpmerged from fastpmerged_2A
    output:
      file ("${sampleid}.bam") into bamRG
      file ("${sampleid}.stats")
    script:
    """
      bwa mem -M -k 10 -t ${threads} ${reference} ${fastpmerged} | samtools sort -T temp -O bam -o ${sampleid}.long.bam

      bwa mem -M -k 10 -t ${threads} ${reference} ${trimmed1} ${trimmed2} | samtools sort -T temp -O bam -o ${sampleid}.short.bam

      samtools merge ${sampleid}.merged.bam ${sampleid}.long.bam ${sampleid}.short.bam
      samtools sort ${sampleid}.merged.bam -@ ${threads} -o ${sampleid}.bam
      samtools index ${sampleid}.bam
      samtools stats ${sampleid}.bam > ${sampleid}.stats
    """
}

// add read groups

process "add read groups" {
    tag '2B'
    publishDir outDir + "/mapping", mode: 'copy'
    input:
      file bam from bamRG
    output:
      file ("${sampleid}.RG.bam") into bamMD
    script:
    """
    gatk AddOrReplaceReadGroups I=${bam} O=${sampleid}.RG.bam SORT_ORDER=coordinate RGID=RG1 RGLB=DN912945U RGPU=HYJN7DRXY RGPL=ILLUMINA RGSM=${sampleid} RGPI=150 CREATE_INDEX=True
    """
}

// mark duplicate sequences with picard

process "mark duplicates"{
    tag '2C'
    publishDir outDir + "/mapping", mode: 'copy'
    input:
      file bam from bamMD
    output:
      file ("${sampleid}.markdup.bam") into bamBQR, bamABQR
      file ("${sampleid}.markdup.stats")
    script:
    """
      gatk MarkDuplicates -I ${bam} -O ${sampleid}.markdup.bam -M ${sampleid}.metrics.txt
      samtools stats ${sampleid}.markdup.bam > ${sampleid}.markdup.stats
    """
}

// Base Recalibration Table

process 'base recalibration table' {
    tag '2D'
    publishDir outDir + '/mapping', mode: 'copy'
    input:
      file bam from bamBQR
    output:
      file ("${sampleid}.recal.table") into recalTableBQR
    script:
    """
    gatk BaseRecalibrator -I ${bam} -R ${reference} --known-sites ${known_sites} -O ${sampleid}.recal.table
    """
}


// Apply Base Recalibration

process 'applybase recalibration' {
    tag '2E'
    publishDir outDir + '/mapping', mode: 'copy'
    input:
      file recalTable from recalTableBQR
      file bam from bamABQR
    output:
      file ("${sampleid}.BQR.bam") into bamDepth, bamVar, bamNonCov
      file ("${sampleid}.BQR.bai") into bam_index

    script:
    """
    gatk ApplyBQSR -I ${bam} -R ${reference} --bqsr-recal-file ${recalTable} -O ${sampleid}.BQR.bam
    samtools index ${sampleid}.BQR.bam > ${sampleid}.BQR.bam.bai

    """
}

// call variants with GATK

process "variant calling"{
    tag '3A'
    publishDir outDir + "/variants", mode: 'copy'
    input:
      file bam from bamVar
    output:
      file ("${sampleid}.vcf.gz") into vcfGen
      file ("${sampleid}.vcf.gz.tbi")
    script:
    """
      gatk HaplotypeCaller -I ${bam} -O ${sampleid}.vcf.gz -R ${reference} -ERC GVCF -G StandardAnnotation
    """
}

// genotype called variants

process "genotyping variants" {
  tag '3B'
  publishDir outDir + "/variants", mode: "copy"
  input:
    file vcf from vcfGen
  output:
    file ("${sampleid}.gen.vcf.gz") into vcfSelect, vcfNonCov, vcfCon
  script:
  """
  gatk IndexFeatureFile -I ${vcf}

  gatk GenotypeGVCFs -R ${reference} -V ${vcf} -O ${sampleid}.gen.vcf.gz --annotation-group StandardAnnotation
  """

}

// subset snps and indels

process 'subset snps and indels' {
  tag '3C'
  publishDir outDir + "/variants", mode: "copy"
  input:
    file vcf from vcfSelect
  output:
    file ("${sampleid}.snp.vcf.gz") into vcfFiltersnp, vcfVQRsnp, vcfAVQRsnp
    file ("${sampleid}.indel.vcf.gz") into vcfFilterindel, vcfVQRindel, vcfAVQRindel
  script:
  """
    gatk IndexFeatureFile -I ${vcf}
    
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include SNP -O ${sampleid}.snp.vcf.gz
    
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include INDEL -O ${sampleid}.indel.vcf.gz
  """
}

/*
// SNPs score recalibration

process "SNP score recalibration" {
  tag '3C'
  publishDir outDir + "/variants", mode: "copy"
  publishDir outDir + "/variants_rscripts", mode: "copy", pattern: "${sampleid}.plots.R"

  input:
    file vcf from vcfVQRsnp
  output:
    file ("${sampleid}snp.recal.table") into SnpRecalTableVQR
    file ("${sampleid}snp.tranches") into SnpTranchesVQR
  script:
  """
    gatk IndexFeatureFile -I ${vcf}
    
    gatk VariantRecalibrator -R ${reference} -V ${vcf} --resource:hapmap,known=true,training=true,truth=true,prior=15.0 ${known_sites} -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR --trust-all-polymorphic true --max-gaussians 8 -mode SNP -O ${sampleid}.snp.recal.table --tranches-file ${sampleid}.snp.tranches --rscript-file ${sampleid}.snp.plots.R
  """
}

// apply SNPs quality recalibration

process "apply SNP score recalibration" {
  tag '3D'
  publishDir outDir + "/variants", mode: "copy"
  input:
  file vcf from vcfAVQRsnp
  file recalTable from SnpRecalTableVQR
  file tranches from SnpTranchesVQR
  output:
    file ("${sampleid}.snp.VQR.vcf.gz") into vcfAnnotsnp
  script:
  """
    gatk ApplyVQSR -R ${reference} -V ${vcf} --truth-sensitivity-filter-level 99.0 --tranches-file ${tranches}  --recal-file ${recalTable} -mode SNP -O ${sampleid}.snp.VQR.vcf.gz
  """
}
*/


// filter SNPs

process "filter SNPs"{
    tag '3F'
    publishDir outDir + "/variants", enabled: false
    input:
      file vcf from vcfFiltersnp
    output:
      file ("${sampleid}.snp.filtered.vcf") into vcfPassedSnps
    script:
    """
      gatk IndexFeatureFile -I ${vcf}

      gatk VariantFiltration -V ${vcf} -O ${sampleid}.snp.filtered.vcf --filter-expression 'QD < ${QD}' --filter-name "LowDepth"  --filter-expression 'FS > ${FS}' --filter-name "FisherStrandBiasFilter" --filter-expression 'MQ < ${MQ}' --filter-name "RMSMappingQualityFilter" --filter-expression 'SOR > ${SOR}' --filter-name "SORFilter" --filter-expression 'ReadPosRankSumTest < ${ReadPosRankSum}' --filter-name "ReadPosRankSumFilter" --filter-expression 'MQRankSum < ${MQRankSum}' --filter-name "MQRankSumFilter"
    """
}

// select only passed SNPs

process 'passed SNPs' {
  tag '3C'
  publishDir outDir + "/variants", mode: "copy"
  input:
    file vcf from vcfPassedSnps
  output:
      file ("${sampleid}.snp.passed.ren.vcf.gz") into vcfAnnotsnp
      file ("${sampleid}.snp.passed.vcf.gz")
      file ("${sampleid}.snp.filtered.vcf.gz")

  script:
  """
    gatk IndexFeatureFile -I ${vcf}

    gatk SelectVariants -R ${reference} -V ${vcf} --exclude-filtered -O ${sampleid}.snp.passed.vcf


    sed -i 's/Pf3D7_01_v3/1/g;\
      s/Pf3D7_02_v3/2/g; \
      s/Pf3D7_03_v3/3/g; \
      s/Pf3D7_04_v3/4/g; \
      s/Pf3D7_05_v3/5/g; \
      s/Pf3D7_06_v3/6/g; \
      s/Pf3D7_07_v3/7/g; \
      s/Pf3D7_08_v3/8/g; \
      s/Pf3D7_09_v3/9/g; \
      s/Pf3D7_10_v3/10/g; \
      s/Pf3D7_11_v3/11/g; \
      s/Pf3D7_12_v3/12/g; \
      s/Pf3D7_13_v3/13/g; \
      s/Pf3D7_14_v3/14/g' \
      ${sampleid}.snp.passed.vcf
      > ${sampleid}.snp.passed.ren.vcf

      sed -i '/Pf3D7_API_v3/d;\
      /Pf3D7_MIT_v3/d'\
      ${sampleid}.snp.passed.ren.vcf
      
      bgzip ${sampleid}.snp.passed.ren.vcf
      
      bgzip ${sampleid}.snp.passed.vcf

      bgzip ${vcf}
    
  """
}

// annotate SNPs with snpEff

process 'annotation SNPs'{
    tag '5A'
    publishDir outDir + "/annotation/snp", mode: 'copy'
    publishDir outDir + '/QC', mode: 'copy', pattern: "${sampleid}.snp.snpEff.csv"
    input:
      file vcf from vcfAnnotsnp
    output:
      file ("${sampleid}.snp.snpEff.csv") into snpEffstatsMQsnps
      file ("${sampleid}.snp.summary.html")
      file ("${sampleid}.snp.snpEff.genes.txt")
      file ("${sampleid}.snp.annot.table.tsv")
      file ("${sampleid}.snp.snpEff.ann.vcf")
    script:
    """
      snpEff ann -v  -strict -noShiftHgvs -s ${sampleid}.snp.summary.html Plasmodium_falciparum -csvStats ${sampleid}.snp.snpEff.csv ${vcf} > ${sampleid}.snp.snpEff.ann.vcf

      ${vcf2table} ${sampleid}.snp.snpEff.ann.vcf --sample ${sampleid} -ad -e -o ${sampleid}.snp.annot.table.tsv


    """
}



/*
// INDEL quality recalibration

process "INDEL score recalibration" {
  tag '3C'
  publishDir outDir + "/variants", mode: "copy"
  publishDir outDir + "/variants_rscripts", mode: "copy", pattern: "${sampleid}.plots.R"

  input:
    file vcf from vcfVQRindel
  output:
    file ("${sampleid}.indel.recal.table") into IndelRecalTableVQR
    file ("${sampleid}.indel.tranches") into IndelTranchesVQR
  script:
  """
    gatk IndexFeatureFile -I ${vcf}
    
    gatk VariantRecalibrator -R ${reference} -V ${vcf} --resource:hapmap,known=true,training=true,truth=true,prior=15.0 ${known_sites}  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR --trust-all-polymorphic true --max-gaussians 8 -mode INDEL -O ${sampleid}.indel.recal.table --tranches-file ${sampleid}.indel.tranches --rscript-file ${sampleid}.indel.plots.R
  """
}


// apply INDEL quality recalibration

process "apply INDEL score recalibration" {
  tag '3D'
  publishDir outDir + "/variants", mode: "copy"
  input:
  file vcf from vcfAVQRindel
  file recalTable from IndelRecalTableVQR
  file tranches from IndelTranchesVQR
  output:
    file ("${sampleid}.indel.VQR.vcf.gz") into vcfAnnotindel
  script:
  """
    gatk ApplyVQSR -R ${reference} -V ${vcf} --truth-sensitivity-filter-level 99.0 --tranches-file ${tranches}  --recal-file ${recalTable} -mode INDEL -O ${sampleid}.indel.VQR.vcf.gz
  """
}
*/

// filter INDELs

process "filter INDELs"{
    tag '3E2'
    publishDir outDir + "/variants", enabled: false
    input:
      file vcf from vcfFilterindel
    output:
      file ("${sampleid}.indel.filtered.vcf") into vcfPassedIndel
    script:
    """
      gatk IndexFeatureFile -I ${vcf}

    gatk VariantFiltration -V ${vcf} -O ${sampleid}.indel.filtered.vcf --filter-expression 'QD < ${QD}' --filter-name "LowDepth"  --filter-expression 'FS > ${FS}' --filter-name "FisherStrandBiasFilter" --filter-expression 'MQ < ${MQ}' --filter-name "RMSMappingQualityFilter" --filter-expression 'SOR > ${SOR}' --filter-name "SORFilter" --filter-expression 'ReadPosRankSum < ${ReadPosRankSum}' --filter-name "ReadPosRankSumTestFilter" --filter-expression 'MQRankSum < ${MQRankSum}' --filter-name "MQRankSumFilter"
    """
}

// select only passed INDELs

process 'passed INDELs' {
  tag '3C'
  publishDir outDir + "/variants", mode: "copy"
  input:
    file vcf from vcfPassedIndel
  output:
      file ("${sampleid}.indel.passed.ren.gz") into vcfAnnotIndel
      file ("${sampleid}.indel.passed.vcf.gz")
      file ("${sampleid}.indel.filtered.vcf.gz")
  script:
  """
    gatk IndexFeatureFile -I ${vcf}

    gatk SelectVariants -R ${reference} -V ${vcf} --exclude-filtered -O ${sampleid}.indel.passed.vcf


    sed -i 's/Pf3D7_01_v3/1/g;\
      s/Pf3D7_02_v3/2/g; \
      s/Pf3D7_03_v3/3/g; \
      s/Pf3D7_04_v3/4/g; \
      s/Pf3D7_05_v3/5/g; \
      s/Pf3D7_06_v3/6/g; \
      s/Pf3D7_07_v3/7/g; \
      s/Pf3D7_08_v3/8/g; \
      s/Pf3D7_09_v3/9/g; \
      s/Pf3D7_10_v3/10/g; \
      s/Pf3D7_11_v3/11/g; \
      s/Pf3D7_12_v3/12/g; \
      s/Pf3D7_13_v3/13/g; \
      s/Pf3D7_14_v3/14/g' \
      ${sampleid}.indel.passed.vcf
      > ${sampleid}.indel.passed.ren.vcf

      sed -i '/Pf3D7_API_v3/d;\
      /Pf3D7_MIT_v3/d'\
      ${sampleid}.snpindel.passed.ren.vcf
      
      bgzip ${sampleid}.indel.passed.ren.vcf
      
      bgzip ${sampleid}.indel.passed.vcf

      bgzip ${vcf}
  """
}

// annotate INDELs snpEff

process 'annotation INDELs'{
    tag '5B'
    publishDir outDir + "/annotation/indel", mode: 'copy'
    publishDir outDir + '/QC', mode: 'copy', pattern: "${sampleid}.indel.snpEff.csv"
    input:
      file vcf from vcfAnnotIndel
    output:
      file ("${sampleid}.indel.snpEff.csv") into snpEffstatsMQindel
      file ("${sampleid}.indel.summary.html")
      file ("${sampleid}.indel.snpEff.genes.txt")
      file ("${sampleid}.indel.annot.table.tsv")
      file ("${sampleid}.indel.snpEff.ann.vcf")
    script:
    """
      snpEff ann -v -htmlStats ${sampleid}.indel.summary.html Plasmodium_falciparum -csvStats ${sampleid}.indel.snpEff.csv ${vcf} > ${sampleid}.indel.snpEff.ann.vcf

      ${vcf2table} ${sampleid}.indel.snpEff.ann.vcf --sample ${sampleid} -ad -e -o ${sampleid}.indel.annot.table.tsv
      
    """
}

// non-covered regions

process 'non covered regions'{
    tag '4'
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
      file bam from bamNonCov
      file vcf from vcfNonCov
    output:
      file("${sampleid}_noncov.bed")
      file("${sampleid}.uncovered.vcf.gz") into vcfCon2
    script:
    """
      bedtools genomecov -ibam ${bam} -bga | \
      awk '\$4 < ${DP}' | \
      awk '{{print(\$1 \"\\t\" \$2 + 1 \"\\t\" \$3 \"\\tlow_coverage\")}}' |\
      bedtools subtract -a - -b ${vcf} > ${sampleid}_noncov.bed
      
      gatk IndexFeatureFile -I ${vcf}

      if [ ! -s "${sampleid}_noncov.bed" ]; then
        continue

      else 
        bcftools view -R ${sampleid}_noncov.bed -O z -o ${sampleid}.uncovered.vcf.gz ${vcf}
      fi  

    """
}

// consensus calling

process 'consensus'{
    tag '6'
    publishDir outDir + "/consensus", mode: 'copy'
    input:
      file vcf from vcfCon
      file vcf2 from vcfCon2
    output:
      file("${sampleid}.consensus.fasta") into consensus
    script:
        """
        gatk IndexFeatureFile -I ${vcf}

      if [ -e ${vcf2} ]; then

        gatk IndexFeatureFile -I ${vcf2}

        gatk FastaAlternateReferenceMaker -R ${reference} -O ${sampleid}.consensus.fasta -V ${vcf} --snp-mask ${vcf2}
      else
        
        gatk IndexFeatureFile -I ${vcf2}

        gatk FastaAlternateReferenceMaker -R ${reference} -O ${sampleid}.consensus.fasta -V ${vcf}
      fi
  
        sed -i 's/>Pf3D7/>${sampleid}/; s/_v3/''/g' ${sampleid}.consensus.fasta
        """
}

// genome completeness/stats calculation

process 'genome_stats'{
    tag '7'
    publishDir outDir + '/consensus', mode: 'copy'
    input:
        file consensus from consensus
    output:
        file ("${sampleid}.stats.txt")
    script:
        """
        faCount ${consensus} > ${sampleid}.stats.txt
        ${parse_stats} --stats ${sampleid}.stats.txt
        """
}

// Genome depth with mosdepth 

process 'genome depth' {
    tag '8'
    publishDir outDir + '/QC', mode: 'copy'
    input:
      file bam from bamDepth
      file bamindex from bam_index
    output:
      file "*"
    script:
    """
      mosdepth --threads $threads ${sampleid} ${bam}
    """
}

// multiqc

process 'multiQC'{
  tag '9'
  publishDir outDir + '/QC', mode: 'copy'
  input:
    file snpEffStats from snpEffstatsMQsnps
    file snpEffStats from snpEffstatsMQindel
  output:
    file "*"
    file ".command.*"
  script:
  """
    multiqc ${outDir}/QC
  """
}

// run parameters

process 'run parameters'{
    tag '10'
    publishDir outDir, mode: 'copy'
    input:
    output:
    script:
        """
        touch parameters.tsv
        echo "Parameter\tValue" >> parameters.txt
        echo "Coverage"\t "$params.DP"  >> parameters.txt
        echo "Quality By Depth"\t "$params.QD"  >> parameters.txt
        echo "Fisher Strand Bias"\t "$params.FS"  >> parameters.txt
        echo "Median Mapping Quality"\t "$params.MMQ" >> parameters.txt
        echo "RMS Mapping Quality"\t "$params.MQ"  >> parameters.txt 
        """
}