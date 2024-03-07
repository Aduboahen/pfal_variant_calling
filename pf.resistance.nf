#!/usr/bin/env nextflow

nextflow.enable.dsl=1

// parameters to pass to command lineage

params.reads       = "$baseDir/test_fq/test_{1,2}.fastq.gz"
reference          = "$baseDir/reference/Pf3D7.fasta"
sampleid           = "$params.sampleid"
vcf2table          = "$baseDir/scripts/vcf2table.py"
parse_stats        = "$baseDir/scripts/parse_stats.py"
outDir             = "$params.outDir"
params.threads     = 50 
threads            = "$params.threads"
params.varqual	   = 15
varqual		   = "$params.varqual" 
params.varthres    = 0.80
varthres           = "$params.varthres"
params.mapqual     = 15 
mapqual            = "$params.mapqual" 
params.depth       = 5
depth              = "$params.depth"


// channel to get reads as tuples

Channel
    .fromFilePairs(params.reads)
    .ifEmpty{"no such files"}
    .set{reads_ch}

log.info"""

codeBase  : "$baseDir"
sample    : "$params.sampleid"
reads     : "$params.reads"
outDir    : "$params.outDir"
reference : "$baseDir/reference/Pf3D7.fasta"
threads   : "$params.threads"
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
      file ("${sampleid}.bam")  into bam1A
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

// mark duplicate sequences with picard

process "mark duplicates"{
    tag '2B'
    publishDir outDir + "/mapping", mode: 'copy'
    input:
      file bam from bam1A
    output:
      file ("${sampleid}.markdup.bam") into (bamDepth, bamVar)//, bamNonCov)
      file ("${sampleid}.markdup.stats")
      file ("${sampleid}.markdup.bam.bai") into bam_index
    script:
    """
      picard MarkDuplicates -I ${bam} -O ${sampleid}.markdup.bam -M ${sampleid}.metrics.txt
      samtools index ${sampleid}.markdup.bam > ${sampleid}.markdup.bam.bai
      samtools stats ${sampleid}.markdup.bam > ${sampleid}.markdup.stats
    """
}

// Process 2C: Genome depth
process 'genome depth' {
    tag '2C'
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

// call variants with bcftools

process "variant calling"{
    tag '3A'
    publishDir outDir + "/variants", mode: 'copy'
    input:
      file bam from bamVar
    output:
      file ("${sampleid}.vcf") into vcfFilter
    script:
    """
      bcftools mpileup -a DP -B -O u -m 4 -f \
      ${reference} --threads ${threads} ${bam}\
      | bcftools call -mv -O v -o \
      ${sampleid}.vcf --threads ${threads} 
    """
}
 
process"filter variants"{
    tag '3B'
    publishDir outDir + "/variants", mode: 'copy'
    input:
      file vcf from vcfFilter
    output:
      //file ("${sampleid}.filtered.vcf") into (vcfCon)//,vcfNonCov
      file ("${sampleid}.filtered_DP.vcf") into vcfCon
      file ("${sampleid}.renamed.vcf") into vcfAnnot
    script:
    """
      bcftools filter -i 'type="snp"\
      && QUAL>=${varqual} \
      && FORMAT/DP>${depth} \
      && MQ>=${mapqual} \
      && DP4[2]/(DP4[2]+DP4[0])>=${varthres}\
      && DP4[3]/(DP4[3]+DP4[1])>=${varthres}'\
      -g10 -G10 \
      --threads ${threads} \
      -o ${sampleid}.filtered_DP.vcf ${sampleid}.vcf  

      sed 's/Pf3D7_01_v3/1/g;\
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
      ${sampleid}.filtered_DP.vcf \
      > ${sampleid}.renamed.vcf

      sed -i '/Pf3D7_API_v3/d;\
      /Pf3D7_MIT_v3/d'\
      ${sampleid}.renamed.vcf
    """
}
/*
bcftools filter -i 'type="snp" \
&& QUAL>=30 \
&& MQ>=30 \
&& DP4[2]/(DP4[2]+DP4[0])>=0.80\
&& DP4[3]/(DP4[3]+DP4[1])>=0.80'\
-g10 -G10\
${sampleid}.vcf \
-o ${sampleid}.filtered.vcf



process 'non covered regions'{
    tag '3C'
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
      file bam from bamNonCov
      file vcf from vcfNonCov
    output:
      file("${sampleid}_noncov.bed") into noncov
    script:
    """
      bedtools genomecov -ibam ${bam} -bga | \
      awk '\$4 < 5' | \
      awk '{{print(\$1 \"\\t\" \$2 + 1 \"\\t\" \$3 \"\\tlow_coverage\")}}' |\
      bedtools subtract -a - -b ${vcf} > ${sampleid}_noncov.bed
    """
}
*/

// annotate with snpEff

process 'annotation'{
    tag '4A'
    publishDir outDir + "/annotation", mode: 'copy'
    publishDir outDir + '/QC', mode: 'copy', pattern: "${sampleid}.snpEff.csv"
    input:
      file vcf from vcfAnnot
    output:
      file ("${sampleid}.snpEff.csv")
      file ("${sampleid}.summary.html")
      file ("${sampleid}.snpEff.genes.txt")
      file ("${sampleid}.annot.table.tsv")
      file ("${sampleid}.snpEff.ann.vcf")
    script:
    """
      snpEff ann -v -noShiftHgvs -ud 0 -strict -hgvs1LetterAa -s \
      ${sampleid}.summary.html Plasmodium_falciparum\
      -csvStats ${sampleid}.snpEff.csv ${vcf} > ${sampleid}.snpEff.ann.vcf

      ${vcf2table} ${sampleid}.snpEff.ann.vcf --sample ${sampleid} -ad -e -o\
      ${sampleid}.annot.table.tsv
    """
}

// Consensus calling
process "consensus"{
    tag '5A'
    publishDir outDir + "/consensus", mode: 'copy'
    input:
    	file vcf from vcfCon
      //file noncov from noncov
    output:
    	file("${sampleid}.consensus.fasta") into consensus
    script:
      """
    	bcftools view -O z -o ${sampleid}.vcf.gz ${vcf}

    	bcftools index ${sampleid}.vcf.gz

    	bcftools consensus -f ${reference}\
     	 --mark-snv lc\
      	 --mark-del -\
     	 ${sampleid}.vcf.gz  > ${sampleid}.consensus.fasta

    	sed -i 's/>Pf3D7/>${sampleid}/;s/_v3/''/g' ${sampleid}.consensus.fasta
    	"""
}

//--mask ${noncov}\

  // 7B: genome completeness/stats calculation
process 'genome_stats' {
    tag '5B'
    publishDir outDir + '/consensus', mode: 'copy'
    input:
        file consensus from consensus
    output:
        file ("${sampleid}.stats.txt")
    script:
        """
        faCount ${consensus} > ${sampleid}.stats.txt
        $parse_stats --stats ${sampleid}.stats.txt
        """
}
