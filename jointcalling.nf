#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * Developed by the Genome Institute of Singapore for
 * SG10K health / the National Precision Medicine Program Singapore
 *
 * Copyright: 2018 Genome Institute of Singapore
 * License: The MIT License (MIT)
 *
 * See LICENSE for more copyright information
 */

/* Generate a MD5 hash (mariogarcia's solution from
 * https://gist.github.com/ikarius/299062/85b6540c99878f50f082aaee236ef15fc78e527c)
 */
import java.security.MessageDigest
def generateMD5_A(String s){
    MessageDigest.getInstance("MD5").digest(s.bytes).encodeHex().toString()
}

workflow_name = "SG10K health: GRCh38 joint-discovery-gatk4"
log.info "======================================"
log.info " ${workflow_name}"
log.info "======================================"


def helpMessage() {
    log.info """
    Usage: nextflow main.nf -params-file sample.yaml --publishdir outdir -profile nscc
    Options:
    -params-file    Sample config file
    -profile        Config for jobs (use 'nscc' for NSCC Aspire 1)
    --publishDir    Copies the process output files to a specified folder
    --keep_workdir  Don't delete workdir
    """.stripIndent()

}
if (params.help) {
    helpMessage()
    exit 0
}

if (params.keep_workdir) {
   log.warn "Not cleaning up work automatically"
   cleanup = false
}
genome = file(params.references.genome)
genome_index = file(params.references.genome + ".fai")
genome_dict = file(params.references.genome.toString().replace('fasta', 'dict'))
dbsnp = file(params.references.dbsnp)
dbsnp_index = file(params.references.dbsnp) + ".tbi"
hapmap = file(params.references.hapmap)
hapmap_idx = file(params.references.hapmap) + ".tbi"
omni = file(params.references.omni)
omni_idx = file(params.references.omni) + ".tbi"
mills = file(params.references.mills)
mills_index = file(params.references.mills) + ".tbi"
phase1_snps = file(params.references.phase1_snps)
phase1_snps_idx = file(params.references.phase1_snps) + ".tbi"
golden_indel = file(params.references.golden_indel)
golden_indel_idx = file(params.references.golden_indel)+ ".tbi"
calling_interval_list = file(params.references.calling_interval_list)
sample_key = "joint"// output file name

if (params.publishdir == null)
    exit 1, "Missing publishdir param"

assert params.references.genome  != null: "Missing reference genome param"
assert params.sample_name_map != null: "Missing sample-to-gvcf map parameter"
assert dbsnp.exists() != null: "Missing $dbsnp param"
assert hapmap.exists() != null: "Missing $hapmap param"
assert omni.exists() != null: "Missing $omni param"
assert mills.exists() != null: "Missing $mills param"
assert phase1_snps.exists() != null: "Missing $phase1_snps param"
assert golden_indel.exists()!= null: "Missing $golden_indel param"
assert calling_interval_list.exists()!= null: "Missing $calling_interval_list param"

// region channel: key, list
//idx = 1
//Channel.from( params.references.region_clusters  )
//    .map { reg_list -> tuple("r" + idx++, reg_list) }
//    .set { region_list_ch }
//region_list_ch.subscribe { println "region_list_ch $it" }

idx = 1
region_list_ch = Channel
     .fromPath(calling_interval_list)
     .splitText().filter{ ! it.startsWith("@") }
     //.map{ it.split('\t')[0..2].join('\t') }
     .map{ reg_list -> tuple("r" + idx++, reg_list.split('\t')[0..2]) }
     //.subscribe { println it }

// vcf.gz and tbi. without the latter you get ImportGVCFs error
// message "Failed to create reader from file" (I think, AW)
Channel.from(params.sample_name_map.values())
    .flatMap { it -> tuple(file(it), file(it + ".tbi")) }
    .set { gvcfs_ch } 
//gvcfs_ch.subscribe { println "gvcfs_ch $it" }


region_and_gvcf_ch = region_list_ch.combine(gvcfs_ch.collect().toList() )
    //.subscribe { println "region_and_gvcf_ch $it" }


/* final for paraoia */
final sample_name_map_str = params.sample_name_map.
    collect { k,v -> k + "\t" + v.tokenize('/').last() }.join('\n')
/* this get sometimes completely mangled with duplicate rows, mismatching entries etc.
def gen_sample_map_str() {
  str = ""
  params.sample_name_map.each{ sample, path ->
    basename = path.tokenize('/').last();
    str += "${sample}\t${basename}\n" 
  }
  return str
}
*/
/* creating a file instead seems to happen in the wrong directory
  def sample_name_map_file = new File('sample_map.csv')
  sample_name_map_file.withWriter('UTF-8') { writer ->
    params.sample_name_map.each{ sample, path ->
      basename = path.tokenize('/').last();
      writer.write"${sample}\t${basename}" }
  }
*/


// FIXME ideally GenomicsDB would receive split files (pre merge)
// FIXME ideally GenomicsDB uses same regions as main.nf. was there a limitation that only one chrom/bed region works in Genomics DB?
process GenomicsDB {
   input:
   file(genome)
   file(genome_index)
   file(genome_dict)
   set region_no, region_list, file(gvcfs) from region_and_gvcf_ch

   output: 
   // main output is a directory
   set region_no, file("ws-reg-${region_no}") into genomicsdb_ch

   script:
   //sample_name_map_str = gen_sample_map_str()
   bed_str = region_list.join("\t").replace(",", "\t").replace("-", "\t")
   //region_no = generateMD5_A(region_list.toString())[0..8]
   """
   echo "${bed_str}" > reg-${region_no}.bed;
   # printf to avoid newline at end, which GenomicsDBImport is allergic to
   # printf '%s'
   echo "${sample_name_map_str}" > sample_map.csv
   ls -l
   gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:ConcGCThreads=${task.cpus} -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} " \
       GenomicsDBImport --sample-name-map sample_map.csv --genomicsdb-workspace-path ws-reg-${region_no} --intervals reg-${region_no}.bed
   """
}


process GenotypeGVF {
    tag "Running GenotypeGVF on region_no $region_no"

    input:
    set region_no, file("ws-reg-${region_no}") from genomicsdb_ch
    file(genome)
    file(genome_index)
    file(genome_dict)

    output:
    set region_no, file("reg-${region_no}.vcf.gz"), file("reg-${region_no}.vcf.gz.tbi") into region_gtvcf_ch

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:ConcGCThreads=${task.cpus} -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} " \
        GenotypeGVCFs -R ${genome} -V gendb://ws-reg-${region_no} -G StandardAnnotation -O reg-${region_no}.vcf.gz
    """
}


process gtvcf_merge {
    tag "Running gtvcf_merge"
    publishDir "${params.publishdir}", mode: 'copy'

    input:
    //set sample_key, file(gt_vcfs), file(gt_tbis) from region_gtvcf_ch.groupTuple()
    file '*' from region_gtvcf_ch.collect()

    output:
    set sample_key, file("${sample_key}.raw.vcf.gz"), file("${sample_key}.raw.vcf.gz.tbi") into gtvcf_merge_ch1

    script:
    """
    bcftools concat -a -o ${sample_key}.tmp.raw.vcf.gz -O z --threads ${task.cpus} reg*.vcf.gz;
    bcftools index -t ${sample_key}.tmp.raw.vcf.gz;
    picard -Dsamjdk.compression_level=2 -Xms4000m -Xmx${task.memory.toGiga()}G \
        -XX:ConcGCThreads=${task.cpus} -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} \
        SortVcf \
        TMP_DIR=\$(dirname ${sample_key})/tmp \
        I=${sample_key}.tmp.raw.vcf.gz \
        O=${sample_key}.raw.vcf.gz;
    """
}

//vqsr steps
// https://software.broadinstitute.org/gatk/documentation/article.php?id=1259
// VQSR params from https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/1.0.3/joint-discovery-gatk4.wdl
process HardFilterAndMakeSitesOnlyVcf {
    tag "Running HardFilterAndMakeSitesOnlyVcf"
    input:
        set sample_key, file("${sample_key}.raw.vcf.gz"), file("${sample_key}.raw.vcf.gz.tbi") from gtvcf_merge_ch1
    output:
        //set sample_key, file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") into sites_only_vcf_ch1, sites_only_vcf_ch2, sites_only_vcf_ch3
        set sample_key, file("${sample_key}.variant_filtered.vcf.gz"), file("${sample_key}.variant_filtered.vcf.gz.tbi"), file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") into sites_only_vcf_ch1, sites_only_vcf_ch2, sites_only_vcf_ch3
    script:
    """
    gatk VariantFiltration \
        --filter-expression "ExcessHet > 54.69" \
        --filter-name ExcessHet \
        -O ${sample_key}.variant_filtered.vcf.gz \
        -V ${sample_key}.raw.vcf.gz;

    gatk  MakeSitesOnlyVcf \
        --INPUT ${sample_key}.variant_filtered.vcf.gz \
        --OUTPUT ${sample_key}.sites_only_vcf_filename.vcf.gz
    """
}
// max-gaussians reduced from 6 to 4 as sometimes Variant recalibrator fails due to lack of data
//Remove contigs from reference
//https://gatkforums.broadinstitute.org/gatk/discussion/comment/34102
//Remove low depth coverage data. Variant recalibrations as mentioned below
//https://gatkforums.broadinstitute.org/gatk/discussion/3952/variantrecalibrator-no-data-found
// -minNumBad to 1000 
//Issue reported https://gatkforums.broadinstitute.org/gatk/discussion/8450/error-java-lang-illegalargumentexception-no-data-found-when-using-vqsr

process VariantRecalibrator_SNPs {
    tag "sample $sample_key"
    input:
        //set sample_key, file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") from sites_only_vcf_ch1
        set sample_key, file("${sample_key}.variant_filtered.vcf.gz"), file("${sample_key}.variant_filtered.vcf.gz.tbi"), file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") from sites_only_vcf_ch1
        file(genome)
        file(genome_index)
        file(genome_dict)
        file(hapmap)
        file(hapmap_idx)
        file(omni)
        file(omni_idx)
        file(phase1_snps)
        file(phase1_snps_idx)
        file(dbsnp)
        file(dbsnp_index)
    output:
        set sample_key, file("${sample_key}.recalibrate_SNP.recal"), file("${sample_key}.recalibrate_SNP.recal.idx"), file("${sample_key}.recalibrate_SNP.tranches") into variantrecalibrator_SNP_ch
    script:
    """
    gatk VariantRecalibrator \
        -V ${sample_key}.sites_only_vcf_filename.vcf.gz \
        -R $genome \
        --trust-all-polymorphic \
        -resource hapmap,known=false,training=true,truth=true,prior=15.0:./$hapmap \
        -resource omni,known=false,training=true,truth=true,prior=12.0:./$omni \
        -resource 1000G,known=false,training=true,truth=false,prior=10.0:./$phase1_snps \
        -resource dbsnp,known=true,training=false,truth=false,prior=7.0:./$dbsnp \
        -an DP \
        -an QD \
        -an FS \
        -an SOR \
        -an MQ \
        -an MQRankSum \
        -an ReadPosRankSum \
        -mode SNP \
        -an InbreedingCoeff \
        --minimum-bad-variants 1000 \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --max-gaussians 4 \
        --output-model ${sample_key}..snps.model.report \
        -O ${sample_key}.recalibrate_SNP.recal \
        --tranches-file ${sample_key}.recalibrate_SNP.tranches 
    """   
}
//Add -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} 
// for GRCh38 build
process VariantRecalibrator_INDELs {
    tag "sample $sample_key"
    input:
        //set sample_key, file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") from sites_only_vcf_ch2
        set sample_key, file("${sample_key}.variant_filtered.vcf.gz"), file("${sample_key}.variant_filtered.vcf.gz.tbi"), file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi") from sites_only_vcf_ch2
        file(genome)
        file(genome_index)
        file(genome_dict)
        file(dbsnp)
        file(dbsnp_index)
        file(golden_indel)
        file(golden_indel_idx)
    output:
        set sample_key, file("${sample_key}.recalibrate_INDEL.recal"), file("${sample_key}.recalibrate_INDEL.recal.idx"), file("${sample_key}.recalibrate_INDEL.tranches") into variantrecalibrator_INDEL_ch
    script:
    """
    gatk VariantRecalibrator \
        -V ${sample_key}.sites_only_vcf_filename.vcf.gz \
        -R $genome \
        --trust-all-polymorphic \
        --resource mills,known=false,training=true,truth=true,prior=12.0:./$golden_indel \
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:./$dbsnp \
        -an QD \
        -an DP \
        -an FS \
        -an SOR \
        -an MQRankSum \
        -an ReadPosRankSum \
        -mode INDEL \
        -an InbreedingCoeff \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --max-gaussians 4 \
        -O ${sample_key}.recalibrate_INDEL.recal \
        --tranches-file ${sample_key}.recalibrate_INDEL.tranches \
    """   
}

process ApplyVQSR {
    tag "sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}", mode: 'copy'
    input:
        //set sample_key, file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi"), file("${sample_key}.recalibrate_INDEL.recal"),  file("${sample_key}.recalibrate_INDEL.recal.idx"), file("${sample_key}.recalibrate_INDEL.tranches"), file("${sample_key}.recalibrate_SNP.recal"), file("${sample_key}.recalibrate_SNP.recal.idx"), file("${sample_key}.recalibrate_SNP.tranches") from sites_only_vcf_ch3.join(variantrecalibrator_INDEL_ch).join(variantrecalibrator_SNP_ch)
       // SHOULD-LIKEY-BE // set sample_key, file("${sample_key}.variant_filtered.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi"), file("${sample_key}.recalibrate_INDEL.recal"),  file("${sample_key}.recalibrate_INDEL.recal.idx"), file("${sample_key}.recalibrate_INDEL.tranches"), file("${sample_key}.recalibrate_SNP.recal"), file("${sample_key}.recalibrate_SNP.recal.idx"), file("${sample_key}.recalibrate_SNP.tranches") from sites_only_vcf_ch3.join(variantrecalibrator_INDEL_ch).join(variantrecalibrator_SNP_ch)
       set sample_key, file("${sample_key}.variant_filtered.vcf.gz"), file("${sample_key}.variant_filtered.vcf.gz.tbi"), \
            file("${sample_key}.sites_only_vcf_filename.vcf.gz"), file("${sample_key}.sites_only_vcf_filename.vcf.gz.tbi"), \
            file("${sample_key}.recalibrate_INDEL.recal"),  file("${sample_key}.recalibrate_INDEL.recal.idx"), \
            file("${sample_key}.recalibrate_INDEL.tranches"), file("${sample_key}.recalibrate_SNP.recal"), \
            file("${sample_key}.recalibrate_SNP.recal.idx"), file("${sample_key}.recalibrate_SNP.tranches") \
            from sites_only_vcf_ch3.join(variantrecalibrator_INDEL_ch).join(variantrecalibrator_SNP_ch)

        //set sample_key, file("${sample_key}.recalibrate_INDEL.recal"), file("${sample_key}.recalibrate_INDEL.recal.idx"), file("${sample_key}.recalibrate_INDEL.tranches") from variantrecalibrator_INDEL_ch
        //set sample_key, file("${sample_key}.recalibrate_SNP.recal"), file("${sample_key}.recalibrate_SNP.recal.idx"), file("${sample_key}.recalibrate_SNP.tranches") from variantrecalibrator_SNP_ch
    output:
         set sample_key, file("${sample_key}.recalibrated.vcf.gz"), file("${sample_key}.recalibrated.vcf.gz.tbi") into results_ch
    script:
    """
    gatk ApplyVQSR \
      -O ${sample_key}.tmp.indel.recalibrated.vcf.gz \
      -V ${sample_key}.variant_filtered.vcf.gz \
      --recal-file ${sample_key}.recalibrate_INDEL.recal \
      --tranches-file ${sample_key}.recalibrate_INDEL.tranches \
      --truth-sensitivity-filter-level 99.7 \
      --create-output-variant-index true \
      -mode INDEL;

    gatk  ApplyVQSR \
      -O ${sample_key}.recalibrated.vcf.gz \
      -V ${sample_key}.tmp.indel.recalibrated.vcf.gz \
      --recal-file ${sample_key}.recalibrate_SNP.recal \
      --tranches-file ${sample_key}.recalibrate_SNP.tranches \
      --truth-sensitivity-filter-level 99.7 \
      --create-output-variant-index true \
      -mode SNP
    """
}


/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html 
 */
workflow.onComplete {
    def msg = """\
    Pipeline execution summary
    ---------------------------
    Status:      : ${ workflow.success ? 'OK' : 'FAILED' }

    Started at   : ${workflow.start}
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}

    Work Dir     : ${workflow.workDir}
    Launch Dir   : ${workflow.launchDir}
    Project Dir  : ${workflow.projectDir}
    """.stripIndent()

    if (! workflow.success) {
       def errmsg = """\

       Report for task that caused the workflow execution to fail:
       Exit status  : ${workflow.exitStatus}
       Error message : ${workflow.errorMessage}
       Error report : ${ workflow.errorReport ? workflow.errorReport : '-' }
        """.stripIndent()

       msg = msg + errmsg
    } 
    status = workflow.success ? 'completed' : 'failed'
    sendMail(from: 'rpd@gis.a-star.edu.sg', to: "${params.mail_to}", 
             subject: "Nextflow execution ${status}: ${workflow_name}", body: msg)
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
