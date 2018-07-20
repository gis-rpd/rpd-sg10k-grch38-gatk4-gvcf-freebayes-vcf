#!/usr/bin/env nextflow
/*
 * vim: syntax=groovy
 * -*- mode: groovy;-*-
 *
 */

/* How to generate a MD5 hash in Groovy: mariogarcia's solution from
 * https://gist.github.com/ikarius/299062/85b6540c99878f50f082aaee236ef15fc78e527c
 */
import java.security.MessageDigest
def generateMD5_A(String s){
    MessageDigest.getInstance("MD5").digest(s.bytes).encodeHex().toString()
}

workflow_name = "SG10K: GRCh38 GATK4-gVCF Freebayes-VCF"
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
    --gvcf_only     Don't compute (Freebayes) vcf
    --keep_workdir  Don't delete workdir
    """.stripIndent()

}
if (params.help) {
    helpMessage()
    exit 0
}


// Check that Nextflow version is up to date enough
try {
    if ( ! nextflow.version.matches(">= $params.nf_required_version") ) {
        throw GroovyException("Nextflow version too old: ${workflow.nextflow.version} < $params.nf_required_version")
    }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}


/* File setup and parameter checks
 * ----------------------------------------------------------------------
 */
if (params.publishdir == null)
    exit 1, "Missing publishdir param"
if (params.samples == null)
    exit 1, "No samples given"
log.info "List of samples: " +  params.samples.keySet()

ref = file( params.references.genome )
if ( ! ref.exists())
    exit 1, "Missing genome reference file: ${ref}"
ref_fai = file( params.references.genome + '.fai')
ref_dict = file( params.references.genome.toString().replace('fasta', 'dict'))
ref_sa = file( params.references.genome + '.sa' )
ref_bwt = file( params.references.genome + '.bwt' )
ref_ann = file( params.references.genome + '.ann' )
ref_amb = file( params.references.genome + '.amb' )
ref_pac = file( params.references.genome + '.pac' )
dbsnp = file( params.references.dbsnp )
if ( ! dbsnp.exists())
    exit 1, "Missing DBSnp reference file: ${dbsnp}"
dbsnp_index = file( params.references.dbsnp ) + '.tbi'

mills = file( params.references.mills )
if ( ! mills.exists())
    exit 1, "Missing Mills reference file: ${mills}"
mills_index = file( params.references.mills ) + '.tbi'

calling_interval_list = file( params.references.calling_interval_list )
if ( ! calling_interval_list.exists())
    exit 1, "Missing Mills reference file: ${calling_interval_list}"

// FIXME we should really test all of them...


/* Channel setup
 * ----------------------------------------------------------------------
 */

def GetReadPair = { sk, rk ->
    // FIXME if files don't exist, their path might be relative to the input yaml
    // see https://gist.github.com/ysb33r/5804364
    tuple(file(params.samples[sk].readunits[rk]['fq1']),
          file(params.samples[sk].readunits[rk]['fq2']))
}

def GetReadUnitKeys = { sk ->
    params.samples[sk].readunits.keySet()
}

Channel
    .from(params.samples.keySet())
    .flatMap { sk -> GetReadUnitKeys(sk).collect{ tuple(sk, it, GetReadPair(sk, it)).flatten() } }
    .set { readunits_ch }
//readunits_ch.subscribe { println "readunits_ch $it" }


cont_vcf_keys = params.references.cont_vcfs.keySet()
cont_vcfs = Channel
    .from ( cont_vcf_keys )
    .map {
    [it, file(params.references.cont_vcfs[it])]
}
//cont_vcfs.subscribe { log.info "value: $it[0]" }


/* And...go!
 * ----------------------------------------------------------------------
 */

process trimadap_map {
    tag { "Aligning readdunit $ru_key of sample $sample_key" }
    input:
        set sample_key, ru_key, file(fq1), file(fq2) from readunits_ch
        file(ref)
        file(ref_amb)
        file(ref_ann)
        file(ref_bwt)
        file(ref_pac)
        file(ref_sa)
    output:
        // also return sample_key for this readunit, to facilitate later merging
        set sample_key, file("${ru_key}.bam"), file("${ru_key}.bam.bai") into ru_bam_ch
    script:
        // sorting threads affect memory, so be nice
        st_threads = (int) Math.ceil(task.cpus/2)
        sort_mem='250M'
        readgroup = "@RG\\tID:$ru_key\\tPU:$ru_key\\tSM:$sample_key\\tLB:$sample_key\\tPL:ILLUMINA"
        """
        seqtk mergepe ${fq1} ${fq2} | \
         trimadap-mt -p ${task.cpus} - | \
            bwa mem -R \"${readgroup}\" -t ${task.cpus} -p ${ref} - | \
            samtools fixmate -@ ${task.cpus} - - | \
            samtools sort -@ ${st_threads} -m ${sort_mem} -T ${ru_key}.tmp -o ${ru_key}.bam;
        samtools index -@ ${task.cpus} ${ru_key}.bam
        """
}

// merge bams is a waste if we only have one pair as input, but in SG10K we should have min 8
process merge_bams {
    tag "Merging readunit BAMs for sample $sample_key"
    input:
        set sample_key, file(bam), file(bai) from ru_bam_ch.groupTuple()
    output:
        set sample_key, file("${sample_key}.merged.bam") into merged_bam_ch
    script:
        """
        samtools merge -@ ${task.cpus} ${sample_key}.merged.bam ${bam.join(' ')}
        """
}

process mark_duplicates {
    tag "Marking duplicates for sample $sample_key"
    input:
        set sample_key, file("${sample_key}.merged.bam") from merged_bam_ch
    output:
        set sample_key, file("${sample_key}.dedup.bam"),file("${sample_key}.dedup.bam.bai") \
            into dedup_bam1_ch, dedup_bam2_ch
    script:
        // READ_NAME_REGEX usage.. Default value <optimized capture of last three ':'
        // separated fields as numeric values> works for both old and new libraries
        // https://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
        // MarkDuplicates has no multithreading options
        """
        #java -Dsamjdk.compression_level=2 -Xms4000m -jar \${PICARD_JAR}
        picard -Dsamjdk.compression_level=2 -Xms4000m -Xmx${task.memory.toGiga()}G \
            -XX:ConcGCThreads=${task.cpus} -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} \
            MarkDuplicates \
            INPUT=${sample_key}.merged.bam \
            OUTPUT=${sample_key}.dedup.bam \
            METRICS_FILE=${sample_key}.metrics.txt \
            VALIDATION_STRINGENCY=SILENT \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            ASSUME_SORTED="true" \
            CLEAR_DT="false" \
            ADD_PG_TAG_TO_READS=false \
            TMP_DIR=\$(dirname ${sample_key})/tmp;
        samtools index ${sample_key}.dedup.bam
        """
}

process gatk_recalibrate_info {
    tag "Running BaseRecalibrator for sample $sample_key"
    input:
        set sample_key, file("${sample_key}.dedup.bam"), file("${sample_key}.dedup.bam.bai") from dedup_bam1_ch
        file(ref)
        file(ref_fai)
        file(ref_dict)
        file(dbsnp)
        file(dbsnp_index)
        file(mills)
        file(mills_index)
    output:
        set sample_key, file("${sample_key}.bqsr.grp") into bqsr_ch
        // FIXME wdl implementation has scatter gather approach!
    script:
        """
        gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms4000m -XX:ConcGCThreads=${task.cpus} \
            -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus}" \
            BaseRecalibrator -R ${ref} -I ${sample_key}.dedup.bam \
            --known-sites ${dbsnp} --known-sites ${mills} -O ${sample_key}.bqsr.grp
        """
}

process gatk_recalibrate_bam {
    tag "Running ApplyBQSR for sample $sample_key"
    input:
        set sample_key, file("${sample_key}.dedup.bam"), file("${sample_key}.dedup.bam.bai"), file("${sample_key}.bqsr.grp")\
            from dedup_bam2_ch.join(bqsr_ch)
        file(ref)
        file(ref_fai)
        file(ref_dict)
    output:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") into \
            bqsr_bam_ch1, bqsr_bam_ch2, bqsr_bam_ch3, bqsr_bam_ch4, bqsr_bam_ch5
        file("${sample_key}.bqsr.bam") into bqsr_bam_only_ch
        file("${sample_key}.bqsr.bam.bai") into  bqsr_bai_only_ch
    script:
        """
        # keeping original qualities, so that we can always go back
        gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:ConcGCThreads=${task.cpus} \
            -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} -Dsamjdk.compression_level=2" \
            ApplyBQSR -R ${ref} -I ${sample_key}.dedup.bam --emit-original-quals \
            -bqsr ${sample_key}.bqsr.grp -O ${sample_key}.bqsr.bam;
        samtools index -@ ${task.cpus} ${sample_key}.bqsr.bam
        """
}

process indexcov {
    tag "Running indexcov for all samples"
    publishDir "${params.publishdir}", mode: 'copy'
    input:
        file bams from bqsr_bam_only_ch.collect()
        file bais from bqsr_bai_only_ch.collect()
    output:
        file("indexcov/all*") into indexcov_ch
    script:
        """
        goleft indexcov -d indexcov/all ${bams}
        """
}

process sample_qc {
    tag "Running Stats/QC for sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}", mode: 'copy'
    input:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") from bqsr_bam_ch1
    output:
        set sample_key, file("${sample_key}.cov-062017.txt"), file("stats/*") into sample_qc_ch
    script:
        """
        sg10k-cov-062017.sh ${sample_key}.bqsr.bam > ${sample_key}.cov-062017.txt;
        mkdir stats;
        samtools stats ${sample_key}.bqsr.bam > stats/${sample_key}.stats;
        plot-bamstats -p stats/${sample_key} stats/${sample_key}.stats
        """
}

process verifybamid {
    tag "Running contamination checks for sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}/verifybamid/", mode: 'copy'
    input:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") from bqsr_bam_ch2
        each file(cont_vcfs)
    output:
        set sample_key, file("*") into verifybamid_ch
    script:
        cont_vcfs = cont_vcfs[1]
        out_base = cont_vcfs.toString().replace(".sites.vcf.GRCh38.liftover.gz", " ")
    """
    verifyBamID --vcf ${cont_vcfs} --bam ${sample_key}.bqsr.bam --out ${sample_key}.${out_base} \
        --noPhoneHome --precise --maxDepth 100 --minMapQ 20 \
        --minQ 20 --maxQ 100
    """
}

process bam2cram {
    tag "Converting BAM to CRAM for sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}", mode: 'copy'
    input:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") from bqsr_bam_ch4
        file(ref)
    output:
        set sample_key, file("${sample_key}.bqsr.cram"), file("${sample_key}.bqsr.cram.crai")
    script:
        """
        samtools view -@ ${task.cpus} -C ${sample_key}.bqsr.bam -o ${sample_key}.bqsr.cram -O CRAM -T ${ref};
        samtools index ${sample_key}.bqsr.cram
        """
}

process freebayes {
    tag "Calling variants with Freebayes for sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}", mode: 'copy'
    input:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") from bqsr_bam_ch5
        file(ref)
        file(ref_fai)
        file(calling_interval_list)
    output:
        set sample_key, file("${sample_key}.fb.vcf.gz"), file("${sample_key}.fb.vcf.gz.tbi")
    when:
        ! params.gvcf_only
    script:
        """
        awk '/^[^@]/ {printf "%s:%d-%d\\n", \$1, \$2, \$3}' ${calling_interval_list} > calling.regions;
        freebayes-parallel calling.regions ${task.cpus} -f ${ref} ${sample_key}.bqsr.bam | bgzip > ${sample_key}.fb-raw.vcf.gz;
        bcftools view -e 'Q<20' -O z -o ${sample_key}.fb.vcf.gz ${sample_key}.fb-raw.vcf.gz;
        tabix -p vcf ${sample_key}.fb.vcf.gz
        """
}

process gatk_hc {
    tag "Running HaplotypeCaller on region $region_no for sample $sample_key"
    //contamination is set to 0 and --max-alternate-alleles set to 6 as f=default (as per documentation)
    //Java option set as per https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.json
    input:
        set sample_key, file("${sample_key}.bqsr.bam"), file("${sample_key}.bqsr.bam.bai") from bqsr_bam_ch3
        each region_list from params['references']['region_clusters']
        file(ref)
        file(ref_fai)
        file(ref_dict)
        file(dbsnp)
        file(dbsnp_index)
    output:
        set sample_key, region_no, file("reg-${region_no}.bed"), \
        file("reg-${region_no}.g.vcf.gz"), file("reg-${region_no}.g.vcf.gz.tbi") into region_gvcf_ch1
    script:
        bed_str = region_list.join("\n").replace(":", "\t").replace("-", "\t")
        region_no = generateMD5_A(region_list.toString())[0..8]
        """
        echo "${bed_str}" > reg-${region_no}.bed;
        gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:ConcGCThreads=${task.cpus}  \
            -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus}" HaplotypeCaller \
            -contamination 0 --max-alternate-alleles 6 -R ${ref} --dbsnp ${dbsnp} \
            -I ${sample_key}.bqsr.bam -L reg-${region_no}.bed --emit-ref-confidence GVCF \
            -O reg-${region_no}.g.vcf.gz
        """
}

process gvcf_merge {
    tag "Merging gVCFs for sample $sample_key"
    publishDir "${params.publishdir}/${sample_key}", mode: 'copy'
    input:
        set sample_key, region_no, file(regbeds), \
            file(reggvcfs), file(reggvcfs_index) from region_gvcf_ch1.groupTuple()
    output:
        file("${sample_key}.g.vcf.gz")
        file("${sample_key}.g.vcf.gz.tbi")
    script:
       """
       bcftools concat -o ${sample_key}.g.vcf.gz -O z --threads ${task.cpus} ${reggvcfs.join(' ')};
       bcftools index --threads ${task.cpus} -t ${sample_key}.g.vcf.gz
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
    Status:      : ${ workflow.success ? 'COMPLETED' : 'FAILED' }

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
       Exit status   : ${workflow.exitStatus}
       Error message : ${workflow.errorMessage}
       Error report  : ${ workflow.errorReport ? workflow.errorReport : '-' }
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
