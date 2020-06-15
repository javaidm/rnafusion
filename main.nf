#!/usr/bin/env nextflow
nextflow.preview.dsl=2

if (params.help) exit 0, helpMessage()

params.name = 'Children Colorado RNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>

// Starting step of the pipeline
stepList = defineStepList()
step = params.step ? params.step.toLowerCase() : ''
if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step ${step}, see --help for more information"

toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'

skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

ch_input_sample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    switch (step) {
        case 'mapping': ch_input_sample = extractFastq(tsvFile); break
        default: exit 1, "Unknown step ${step}"
    }
}

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.fasta_gz = params.genome ? params.genomes[params.genome].fasta_gz ?: null : null
params.fasta_gz_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_gz_fai ?: null : null
params.fasta_gzi = params.genome ? params.genomes[params.genome].fasta_gzi ?: null : null

// The rest can be sorted
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null

// Create the channels for the parmas so Nextflow can stage the files when running the processes
ch_fasta = params.fasta Channel.value(file(params.fasta)) : "null"
ch_intervals = params.intervals && !params.no_intervals  ? Channel.value(file(params.intervals)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"
ch_bait_bed = params.bait_bed ? Channel.value(file(params.bait_bed)) : "null"
ch_fasta_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai))
ch_dict = params.dict ? Channel.value(file(params.dict)) : "null"

printSummary()

/* Check if the fastq needs to be split into multiple files using the Nextflow splitFastQ operator */
// ch_input_pair_reads = Channel.empty()
// if (params.split_fastq){
//         // newly splitfastq are named based on split, so the name is easier to catch
//     ch_input_pair_reads = ch_input_sample
//         // .splitFastq(by: ÷, compress:true, file:"split", flat: true, compress: true, pe:true)
//         .splitFastq(by: 1_000_000, compress:true, file:"split", pe:true)
//         .map {idPatient, idSample, idRun, reads1, reads2 ->
//             // The split fastq read1 is the 4th element (indexed 3) its name is split_3
//             // The split fastq read2's name is split_4
//             // It's followed by which split it's acutally based on the mother fastq file
//             // Index start at 1
//             // Extracting the index to get a new IdRun
//             splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
//             newIdRun = idRun + "_" + splitIndex
//             // Giving the files a new nice name
//             newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
//             newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
//             [idPatient, idSample, newIdRun, reads1, reads2]}
// }
// if (params.split_fastq){
//     ch_input_pair_reads.
//     subscribe{println (it)}
// }



workflow{

    GetSoftwareVersions() 
    BuildIntervals(ch_fasta_fai)
    ch_intervals = params.no_intervals ? "null" : \
                    params.intervals && !('annotate' in step) ? \
                    Channel.value(file(params.intervals)) : BuildIntervals.out
    CreateIntervalBeds(ch_intervals)
    // CreateIntervalBeds.out.flatten().subscribe{println(it)}
    ch_bed_intervals = sortBedIntervalsByDescendingDuration(
                            CreateIntervalBeds.out.flatten() )
    
    if (params.no_intervals  bedIntervals = Channel.from(file("no_intervals.bed"))
    FastQCFQ(ch_input_sample)
   

} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}


/******************************************************************************************/
                                /* Helper functions */
/******************************************************************************************/

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}
// Layer Lab ascii art
   def layerLabAscii() {
    // def ascii_str = 
    return
    '''
    | |                         | |         | |    
    | |     __ _ _   _  ___ _ __| |     __ _| |__ 
    | |    / _` | | | |/ _ \\ '__| |    / _` | '_ \\ 
    | |___| (_| | |_| |  __/ |  | |___| (_| | |_) |
    |______\\__,_|\\__, |\\___|_|  |______\\__,_|_.__/ 
                __/ |                            
               |___/  
    '''
    
  }

def layerLabMessage() {
  // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
     return """ ${c_dim}----------------------------------------------------${c_reset}
     ${c_cyan} Children Colorado RNA Seq ANALYSIS PIPELINE ${c_reset}
     ${c_dim}----------------------------------------------------${c_reset}
     """
  
}

def helpMessage() {
  // Display help message
    log.info layerLabMessage()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --samples samples_dir -profile fiji

    Mandatory arguments:
        --samples                   directory containing samples to be analyzed
        -profile                    Configuration profile to use
                                    Can use multiple (comma separated)
                                    Available: fiji, singularity, test
        --samples                  directory containing samples to be analyzed                                     
        --readlength"                   read length in basepairs                   

    Options:
        --no_intervals              Disable usage of intervals
        --step                      Specify starting step
                                    Available: Mapping
                                    Default: Mapping
        --tools                     Specify tools to be use for analysis
                                    Available: starfusion, arriba, ericscript, squid, suppa
                                    Default: None
        --skip_qc                   Specify which QC tools to skip when running Sarek
                                    Available: all, FastQC, MultiQC, samtools, versions
                                    Default: None
        --target_bed                Target BED file (aka primary/regions/empirical) for targeted  or whole exome sequencing
        --bait_bed                  Bait BED file (aka covered/captured) for targeted or whole exome sequencing (used for GATK CollectHsMetrics)

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
        --fasta                     fasta reference
        --fasta_fai                 reference index
                                    If none provided, will be generated automatically from the fasta reference
         --starindex                directory containing STAR aligner index  
         --salmonindex              directory containing the Salmon aligner index for SUPPA  
                                    If none provided, will be generated automatically if a germlineResource file is provided
        --intervals                 intervals
                                    If none provided, will be generated automatically from the fasta reference
                                    Use --no_intervals to disable automatic generation
        --gtf                       location of annotation gtf                   
                    
        --blacklist                location of blacklist tsv for Arriba                   
        --sf_source                location of STAR-Fusion source code                                   
        --ctat                     directory containing CTAT genome library for STAR-Fusion                                                    
                        
        --fr_db                    directory containing the Fusion Report database                   
        --fusioninspector          perform FusionInspector analysis                   
        --fi_source                location of FusionInspector source code                                   
        --refflat"                  location of refFlat file for (GATK) Picard 
    Other options:
        --outdir                    The output directory where the results will be saved
        --sequencing_center         Name of sequencing center to be displayed in BAM file
        --multiqc_config            Specify a custom config file for MultiQC
        --monochrome_logs           Logs will be without colors
        --email                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --max_multiqc_email_file_size   Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        --exome                     Specify when dealing with WES data, used when calling germline variants using Google deepvariant
        -name                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    AWSBatch options:
        --awsqueue                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion                 The AWS Region for your AWS Batch job to run on
    """.stripIndent()
    // println(params.genome)
}

/*
================================================================================
                                PROCESSES
================================================================================
*/
/*
 * Parse software version numbers
 */
process GetSoftwareVersions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    output:
        // file 'software_versions_mqc.yaml', emit: yamlSoftwareVersion
        file 'software_versions_mqc.yaml'

    when: !('versions' in skipQC)

    script:
    """
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    echo "SNPEFF version"\$(snpEff -h 2>&1) > v_snpeff.txt
    fastqc --version > v_fastqc.txt 2>&1 || true
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    R --version &> v_r.txt  || true
    samtools --version &> v_samtools.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)
process BuildIntervals {
  tag {fastaFai}

  publishDir params.outdir

  input:
    file(fastaFai)

  output:
    file("${fastaFai.baseName}.bed")

  when: !(params.intervals) && !('annotate' in step) && !(params.no_intervals)

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

process CreateIntervalBeds {
    tag {intervals.fileName}

    input:
        file(intervals)

    output:
        file '*.bed'

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotides_per_second}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}


process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip}")

    when: !('fastqc' in skipQC)
    
    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process MapReads_Arriba{
    label 'STAR_Arriba'
    label 'cpus_32'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Mapped/${idSample}/Arriba", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${sample_run}_R1.fastq.gz"), file("${sample_run}_R2.fastq.gz")
        file(genome_dir)
        val(read_length)

    output:
        tuple val('Arriba'), idPatient, idSample, idRun, file("${sample_run}.bam") , emit: bam_mapped
        // tuple val('Arriba'), idPatient, val("${sample_run}"), file("${sample_run}.bam"), emit: bam_mapped_bam_qc
    
    when: 'arriba' in tools
    script:
    sample_run = "${idSample}_${idRun}" 
    overhang = read_length - 1
    """
    STAR \
        --genomeDir $genome_dir \
        --runThreadN ${task.cpus} \
        --readFilesIn  ${sample_run}_R1.fastq.gz ${sample_run}_R2.fastq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 3 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM SoftClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreMin 1 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentReadGapMax 3 \
        --readFilesCommand zcat \
        --sjdbOverhang ${overhang} \
        --outFileNamePrefix ${sample_run}
    """   
}

process MapReads_STAR_Fusion{
    label 'cpus_32'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Mapped/${idSample}/STAR_Fusion", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${sample_run}_R1.fastq.gz"), file("${sample_run}_R2.fastq.gz")
        file(genome_dir)
        val(read_length)

    output:
        tuple val('STAR_Fusion'), idPatient, idSample, idRun, file("${sample_run}.bam") , emit: bam_mapped
        tuple idSample,file("${sample_run}_R1.fastq.gz"), file("${sample_run}_R2.fastq.gz"), file("*_Chimeric.out.junction") , emit: for_star_fusion
        // tuple val('STAR_Fusion'), idPatient, val("${sample_run}"), file("${sample_run}.bam"), emit: bam_mapped_bam_qc
    
    when: 'star_fusion' in tools
    script:
    sample_run = "${idSample}_${idRun}" 
    overhang = read_length - 1
    """
    STAR \
        --genomeDir $genome_dir \
        --runThreadN ${task.cpus} \
        --readFilesIn  ${sample_run}_R1.fastq.gz ${sample_run}_R2.fastq.gz \
        --twopassMode Basic \
        --outReadsUnmapped None \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --chimSegmentReadGapMax 3 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --outSAMattrRGline ID:GRPundef \
        --chimMultimapScoreRange 10 \
        --chimMultimapNmax 10 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --readFilesCommand zcat \
        --chimOutJunctionFormat 1 \
        --sjdbOverhang ${overhang} \
        --outFileNamePrefix ${sample_run}

    """   
}

process MapReads_Squid{
    label 'cpus_32'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Mapped/${idSample}/Squid", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${sample_run}_R1.fastq.gz"), file("${sample_run}_R2.fastq.gz")
        file(genome_dir)
        val(read_length)
        file(gtf)

    output:
        tuple val('Squid'), idPatient, idSample, idRun, file("${sample_run}.bam") , emit: bam_mapped
        tuple idPatient, idSample, idRun, file("${sample_run}.bam"), file("${sample_run}_Chimeric.out.bam") , emit: for_squid
        // tuple val('Squid'), idPatient, val("${sample_run}"), file("${sample_run}.bam"), emit: bam_mapped_bam_qc
    
    when: 'squid' in tools
    script:
    sample_run = "${idSample}_${idRun}" 
    overhang = read_length - 1
    """
    STAR \
        --genomeDir $genome_dir \
        --runThreadN ${task.cpus} \
         --sjdbGTFfile ${gtf} \
        --readFilesIn  ${sample_run}_R1.fastq.gz ${sample_run}_R2.fastq.gz \
        --twopassMode Basic \
        --chimOutType SeparateSAMold \
        --chimSegmentMin 20 \
        --chimJunctionOverhangMin 12 \
        --alignSJDBoverhangMin 10 \
        --outReadsUnmapped Fastx \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --outFileNamePrefix {}/Squid/{}_
        --outFileNamePrefix ${sample_run}

        samtools view -bS ${sample_run}_Chimeric.out.sam > ${sample_run}_Chimeric.out.bam


    """   
}

process IndexBamFile {
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file(bam), file("*.bai")

    script:
    """
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    """
}

process Arriba {
    label 'cpus_16'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Arriba/${idSample}", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file("*.bai")
        file(gtf)
        file(arriba_blacklist)
        fil(fasta)
        fil(fasta_fai)

    output:
        file("${idSample}_arriba.tsv")
        file("${idSample}_discarded_arriba.tsv")

    when: 'arriba' in tools
    
    script:
    """
    arriba \
        -x ${bam} \
        -a ${fasta} \
        -g ${gtf} \
        -b ${arriba_blacklist} \
        -o ${idSample}_arriba.tsv \
        -O ${idSample}_discarded_arriba.tsv \
        -T -P
    """
}

process STAR_Fusion {
    label 'cpus_16'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/STAR_Fusion/${idSample}", 
    mode: params.publish_dir_mode

    input:
        tuple idSample, file(R1), file(R2), file(chimeric_junction)
        file(ctat)

    output:
        file("*.tsv")

    when: 'star_fusion' in tools
    
    script:
    """
   STAR-Fusion \
            --genome_lib_dir ${ctat} \
            -J ${chimeric_junction} \
            --left_fq ${R1} \
            --right_fq ${R2} \
            --CPU ${task.cpus} \
            --examine_coding_effect \ 
            --output_dir STAR_Fusion

            mv STAR_Fusion/star-fusion.fusion_predictions.tsv ${idSample}_star-fusion.tsv
            mv STAR_Fusion/star-fusion.fusion_predictions.abridged.tsv ${idSample}_abridged.tsv
            mv STAR_Fusion/star-fusion.fusion_predictions.abridged.coding_effect.tsv ${idSample}_abridged.coding_effect.tsv
    """
}

process Ericscript {
    label 'cpus_32'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/EricScript/${idSample}", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file(R1), file(R2)
        file(ericscript_db)

    output:
        file("${idSample}_ericscript.tsv")

    when: 'ericscript' in tools
    
    script:
    """
    ericscript.pl \
            -db ${ericscript_db} \
            -name fusions \
            -p ${task.cpus} \
            -o EricScript \
            ${R1} ${R2}
    
    mv EricScript/fusions.results.filtered.tsv ${idSample}_ericscript.tsv
    """
}

process Squid {
    label 'cpus_32'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Squid/${idSample}", 
    mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, idRun, file(bam), file(chimeric_bam)

    output:
        tuple val(idSample), file("${idSample}_fusions_sv.txt")

    when: 'squid' in tools
    
    script:
    """
    squid \
        -b ${bam} \
        -c ${chimeric_bam} \
        -o ${idSample}_fusions_sv.txt
    """
}

process AnnotateSquidOutput {
    label 'cpus_2'
    tag {idSample}

    publishDir "${params.outdir}/Squid_annotated/${idSample}", 
    mode: params.publish_dir_mode

    input:
    tuple val(idSample), file("${idSample}_fusions_sv.txt")
    file(gtf)

    output:
       tuple val(idSample), file("${idSample}_fusions_ann.txt"

    when: 'squid' in tools
    
    script:
    """
    python \
    AnnotateSQUIDOutput.py \
    ${gtf} \
    ${"${idSample}_fusions_sv.txt"} \
    ${idSample}_fusions_ann.txt   
    """
}

process AnnotateSquidOutputPostprocess {
    label 'cpus_2'
    tag {idSample}

    publishDir "${params.outdir}/Postprocessed_Squid_annotations/${idSample}", 
    mode: params.publish_dir_mode

    input:
    tuple val(idSample), file(ann_txt)

    output:
       file("${idSample}_postprocessed_fusions_ann_orig.txt"
       file("${idSample}_postprocessed_fusions_ann.txt"

    when: 'squid' in tools
    
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    
    cols=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', \
                                'name', 'score', 'strand1', 'strand2', 'Type', 'FusedGenes']
    
    fusions = pd.read_csv(${ann_txt}, delimiter='\t', skiprows=1, names = cols)

    if (fusions.FusedGenes.str.split(',').apply(len) > 1).any():
        fusions.to_csv("${idSample}_postprocessed_fusions_annotated_orig.txt", sep='\t', index=False)
        fusions['FusedGenes'] = fusions.FusedGenes.str.split(',').str[0]
        fusions.to_csv("${idSample}_postprocessed_fusions_ann_.txt".format(fastq,fastq.split('/')[-1]),sep='\t',index=False)
    """
}


 
process addEGFRvIII {
    label 'cpus_8'
    tag {idSample}

    publishDir "${params.outdir}/EGFR/${idSample}", 
    mode: params.publish_dir_mode

    input:
    file(gtf)
    file(cdna_fa)
    
    output:
       file("*.cdna.fa")
    
    when 'suppa' in tools

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    """ Modifying gtf file """
    names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    orig = pd.read_csv( $gtf, comment= '#', delimiter= '\t', names= names)
    egfr = orig.loc[ (orig.seqname == 'chr7') & (orig['start'] > 5.5e7) & (orig['end'] < 5.54e7)]
    egfr['gene_id'] = egfr.attribute.str.split('gene_id "').str[1].str.split('";').str[0]
    egfr = egfr.loc[egfr.gene_id.str.contains('ENSG00000146648')]
    egfr['transcript_id'] = egfr.attribute.str.split('; transcript_id "').str[1].str.split('";').str[0]
    egfr['transcript_name'] = egfr.attribute.str.split('; transcript_name "').str[1].str.split('";').str[0]
    egfr['exon_number'] = egfr.attribute.str.split('; exon_number ').str[1].str.split(';').str[0]
    egfr.loc[~egfr.exon_number.isnull(),'exon_number'] = egfr.loc[~egfr.exon_number.isnull(),'exon_number'].astype(int)
    egfr['length'] = egfr['end'] - egfr['start'] + 1
    wt = egfr.loc[egfr.transcript_name == 'EGFR-001']
    if 'EGFRvIII' not in egfr.transcript_id.tolist():
        egfr3 = wt.loc[wt.exon_number.isnull() | (wt.exon_number < 2) | (wt.exon_number > 7)]
        egfr3['attribute'] = 'gene_id "ENSG00000146648.11"; transcript_id "EGFRvIII";'
        egfr3 = egfr3[egfr3.columns[:9]]
        orig = orig.loc[:wt.index.max()].append(egfr3,ignore_index=True)\
        .append(orig.loc[wt.index.max() + 1:],ignore_index=True)
        orig.to_csv(gtf,index=False,header=False,quoting=csv.QUOTE_NONE,sep='\t')
        """ Adding comments """
        tempData = open(gtf,'r')
        orig = tempData.read()
        tempData.close()
        orig = """##description: modified version of GRCh37 to include EGFRvIII
                ##provider: Taylor Firman
                ##contact: taylor.firman@childrenscolorado.org
                ##format: gtf
                ##date: 2020-05-18
                """ + orig
                        tempData = open(gtf,'w')
                        tempData.write(orig)
                        tempData.close()
    """ Modifying cdna file (assuming it has the same filename prefix) """
    tempData = open('.'.join(gtf.split('.')[:-1]) + '.cdna.fa','r')
    sequences = tempData.read()[1:].split('\n>')
    tempData.close()
    transcripts = [seq.split(' ')[0] for seq in sequences]
    if 'EGFRvIII' not in transcripts:
        ind = transcripts.index('ENST00000275493.2')
        egfr = ''.join(sequences[ind].split('\n')[1:])
        egfr3 = egfr[:wt.loc[(wt.feature == 'exon') & (wt.exon_number == 1)].length.sum()] + \
        egfr[wt.loc[(wt.feature == 'exon') & (wt.exon_number < 8)].length.sum():]
        egfr3 = 'EGFRvIII ENSG00000146648.11 EGFR \n' + '\n'.join([egfr3[line*60:(line + 1)*60] for line in range(len(egfr3)//60 + 1)])
        sequences = sequences[:ind + 1] + [egfr3] + sequences[ind + 1:]
        sequences = '>' + '\n>'.join(sequences)
        tempData = open('.'.join(gtf.split('.')[:-1]) + '.cdna.fa','w')
        tempData.write(sequences)
        tempData.close()
    """
}

process Suppa {
    label 'cpus_8'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Suppa/${idSample}", 
    mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, idRun, file("${sample_run}_R1.fastq.gz"), file("${sample_run}_R2.fastq.gz")
    file(salmon_index)
    file(gtf)

    output:
        tuple val(idSample), file("${idSample}_fusions_sv.txt")

    when: 'squid' in tools
    
    script:
    """
    salmon_quant.sh \
            ${salmon_index} \
            "${sample_run}_R1.fastq.gz" \
            "${sample_run}_R2.fastq.gz" \
            ${task.cpus}

    suppa.py generateEvents -i ${gtf} -o ${gtf}.events -e SE -f ioe
    multipleFieldSelection.py -i ${idSample}_quant.sf -k 1 -f 4 -o ${idSample}_iso_tpm.txt
    suppa.py psiPerIsoform -g ${gtf} -e ${idSample}_iso_tpm.txt -o ${idSample}.iso.psi
    """
}

process SuppaPostprocess {
    label 'cpus_8'
    tag {idSample}

    publishDir "${params.outdir}/Suppa/${idSample}", 
    mode: params.publish_dir_mode

    input:
    tuple val(idSample), file(events_psi), file(iso_psi)

    output:
        file("${idSample}.ExonSkipping.txt")

    when: 'suppa' in tools
    
    script:
    """
     #!/usr/bin/env python
    import pandas as pd
    cols = {'index':'event','SUPPA':'PSI'}
    iso = pd.read_csv(${isoform_psi}, delimiter='\t').reset_index().rename(columns=cols),ignore_index=True)
    events = pd.read_csv(${events_psi}, delimiter='\t').reset_index().rename(columns=cols)
    psi = events.append(iso)
    tmp = pd.DataFrame({'event':['ENSG00000146648.11;EGFRvIII',\
    'ENSG00000105976.10;SE:chr7:116411708-116411903:116412043-116414935:+'],\
    'name':['EGFR variant III','MET Exon 14 Skipping']})
    
    psi = pd.merge(left = psi, right = tmp , how = 'inner', on='event')
    psi.to_csv("${idSample}.ExonSkipping.txt")
    """
}

process FusionReport {
    label 'cpus_8'
    tag {idSample}

    publishDir "${params.outdir}/Report/${idSample}", 
    mode: params.publish_dir_mode

    input:
    file(db)
    tuple val(idSample), file(db), file(arriba), file(starfusion), file(squid), file(es)

    output:
        file("Fusion_Report")

    when: 'fusion_report' in tools
    
    script:
    """
    fusion_report run ${idSample} Fusion_Report ${db} \
            --arriba ${arriba} \
            --starfusion ${starfusion} \
            --squid ${squid} \
            --ericscript ${ericscript}
    mv Fusion_Report/fusion_list.tsv Fusion_Report/${idSample}_fusion_list.tsv
    mv Fusion_Report/fusion_genes_mqc.json Fusion_Report/${idSample}_fusion_genes_mqc.json
    """
}

process CollectInsertSizeMetrics{
    label 'cpus_32'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/insert_size_metrics/", mode: params.publish_dir_mode
    
    input:
    tuple val(sampleId), file(bam)

    output:
    file("${sampleId}_aligned_insert_size_metrics.txt")
    file("${sampleId}_aligned_insert_size_histogram.pdf")
    
    when: !('insert_size_metrics' in skipQC)
    script:
    """
    gatk CollectInsertSizeMetrics \
      -I ${bam} \
      -O ${sampleId}_aligned_insert_size_metrics.txt \
      -H ${sampleId}_aligned_insert_size_histogram.pdf \
      M=0.5
    """
}

process CollectRnaSeqMetrics{
    label 'cpus_32'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/insert_size_metrics/", mode: params.publish_dir_mode
    
    input:
    tuple val(sampleId), file(bam)
    file(ref_flat)

    output:
    file("${sampleId}_aligned_rna_metrics.txt")
    
    when: !('rna_seq_metrics' in skipQC)

    script:
    """
    gatk CollectRnaSeqMetrics \
      -I ${bam} \
      -O ${sampleId}_aligned_rna_metrics.txt \
      -REF_FLAT ${ref_flat} \
      -STRAND SECOND_READ_TRANSCRIPTION_STRAND
    """
}

process ResultsExcel{
    label 'cpus_32'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/ResultsExcel/", mode: params.publish_dir_mode
    
    input:
    file(fusion_json)
    file(exon_skipping)

    output:
    file("${sampleId}_detected_fusions.xlsx")
    
    when: !('results_excel' in skipQC)

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    tempData = open(${fusion_json},'r')
    df = pd.DataFrame(json.load(tempData))
    df['Databases'] = df['Databases'].apply(', '.join)
    df["5' Partner Position"] = None
    df["3' Partner Position"] = None
    tempData.close()
    for aligner in ['arriba','starfusion','ericscript','squid']:
        if vars(options)[aligner]:
            inds = ~df[aligner].isnull()
            if inds.sum() > 0:
                for stat in df.loc[inds,aligner].values[0].keys():
                    df.loc[inds,aligner + '_'  + stat] = df.loc[inds,aligner].apply(lambda x: x[stat])
                    if aligner == 'arriba' and stat == 'position':
                        df.loc[inds & df["5' Partner Position"].isnull(),"5' Partner Position"] = \
                        'chr' + df.loc[inds & df["5' Partner Position"].isnull(),aligner + '_' + stat].str.split('#').str[0]
                        df.loc[inds & df["3' Partner Position"].isnull(),"3' Partner Position"] = \
                        'chr' + df.loc[inds & df["3' Partner Position"].isnull(),aligner + '_' + stat].str.split('#').str[1]
                    elif aligner == 'starfusion' and stat == 'position':
                        df.loc[inds & df["5' Partner Position"].isnull(),"5' Partner Position"] = \
                        df.loc[inds & df["5' Partner Position"].isnull(),aligner + '_' + stat].str.split('#').str[0].str[:-2]
                        df.loc[inds & df["3' Partner Position"].isnull(),"3' Partner Position"] = \
                        df.loc[inds & df["3' Partner Position"].isnull(),aligner + '_' + stat].str.split('#').str[1].str[:-2]
            del df[aligner]
    df['Reading Frame'] = df['arriba_reading-frame'] if options.arriba else None
    df['Type'] = df['arriba_type'] if options.arriba else None
    df["5' Partner Coverage"] = df['arriba_coverage1'] if options.arriba else None
    df["3' Partner Coverage"] = df['arriba_coverage2'] if options.arriba else None
    df['Call Confidence'] = df['arriba_confidence'] if options.arriba else None
    df['FFPM'] = df['starfusion_ffmp'] if options.starfusion else None
    df['PSI'] = None
    if options.suppa:
        suppa_vars = pd.read_csv(${exon_skipping}.format(options.samples,options.samples.split('/')[-1]),delimiter='\t')
        suppa_vars = suppa_vars[['name','PSI']].rename(columns={'name':'Fusion'})
        suppa_vars['Type'] = 'exon skipping'
        df = df.append(suppa_vars,ignore_index=True)
    df['Report (Y/N)'] = 'N'
    df = df[['Fusion','Databases','Score','Explained score',"5' Partner Position",\
    "3' Partner Position",'Reading Frame','Type',"5' Partner Coverage",\
    "3' Partner Coverage",'Call Confidence','FFPM','PSI','Report (Y/N)'] + \
    [col for col in df.columns if col.split('_')[0] in ['arriba','starfusion','ericscript','squid']]]
    writer = pd.ExcelWriter("${sampleId}_detected_fusions.xlsx",engine='xlsxwriter')
    f = writer.book.add_format()
    f.set_align('center')
    f.set_align('vcenter')
    df.to_excel(writer,sheet_name='Fusions',index=False)
    for idx, col in enumerate(df):
        series = df[col]
        max_len = min(max((series.astype(str).map(len).max(),len(str(series.name)))) + 1,150)
        writer.sheets['Fusions'].set_column(idx,idx,max_len,f,{'hidden':col.split('_')[0] \
        in ['arriba','starfusion','ericscript','squid'] or col == 'Explained score'})
    writer.sheets['Fusions'].autofilter('A1:' + (chr(64 + (df.shape[1] - 1)//26) + \
    chr(65 + (df.shape[1] - 1)%26)).replace('@','') + str(df.shape[0] + 1))
    writer.sheets['Fusions'].freeze_panes(1,1)
    writer.save()
    """
}


/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

process MultiQC {
    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode
    input:
        file (multiqcConfig) 
        file (versions) 
        // file ('bamQC/*') 
        file ('BCFToolsStats/*') 
        file ('FastQC/*') 
        file ('MarkDuplicates/*') 
        file ('SamToolsStats/*') 
        file ('snpEff/*') 
        file ('VCFTools/*')
        file ('CollectHsMetrics/*')
        file ('CollectAlignmentSummary/*')

    output:
        set file("*multiqc_report.html"), file("*multiqc_data") 

    when: !('multiqc' in skipQC)

    script:
    """
    multiqc -f -v .
    """
}


/******************************************************************************************/
                                /* Helper functions */
/******************************************************************************************/


def printSummary(){
    /*
    ================================================================================
                                    PRINTING SUMMARY
    ================================================================================
    */

    // Header log info
    log.info layerLabMessage()
    def summary = [:]
    if (workflow.revision)          summary['Pipeline Release']    = workflow.revision
    summary['Run Name']          = custom_runName ?: workflow.runName
    summary['Max Resources']     = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
    if (workflow.containerEngine)   summary['Container']         = "${workflow.containerEngine} - ${workflow.container}"
    if (params.samples)               summary['Samples']             = params.input
    if (params.target_bed)           summary['Target BED']        = params.target_bed
    if (params.bait_bed)           summary['BAIT BED']        = params.bait_bed
    if (step)                       summary['Step']              = step
    if (params.tools)               summary['Tools']             = tools.join(', ')
    if (params.skip_qc)              summary['QC tools skip']     = skipQC.join(', ')

    if (params.no_intervals && step != 'annotate') summary['Intervals']         = 'Do not use'
    // if ('strelka' in tools && 'manta' in tools )   summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
    if (params.sequencing_center)                  summary['Sequenced by']      = params.sequencing_center
    if (params.pon && 'mutect2' in tools)          summary['Panel of normals']  = params.pon

    // summary['Save Genome Index'] = params.saveGenomeIndex ? 'Yes' : 'No'
    // summary['Nucleotides/s']     = params.nucleotidesPerSecond
    summary['Output dir']        = params.outdir
    summary['Launch dir']        = workflow.launchDir
    summary['Working dir']       = workflow.workDir
    summary['Script dir']        = workflow.projectDir
    summary['User']              = workflow.userName
    summary['genome']            = params.genome

    if (params.fasta)                 summary['fasta']                 = params.fasta
    if (params.fasta_fai)              summary['fasta_fai']              = params.fasta_fai
    // if (params.fasta_gz)              summary['fasta_gz']              = params.fasta_gz
    // if (params.fasta_gz_fai)          summary['fasta_gz_fai']        = params.fasta_fai
    // if (params.fasta_gzi)              summary['fasta_gzi']              = params.fasta_gzi

    // if (params.dict)                  summary['dict']                  = params.dict
    if (params.intervals)             summary['intervals']             = params.intervals
    if (params.giab_highconf)            summary['GIAB Truth Set']          = params.giab_highconf
    if (params.chco_highqual_snps)     summary['Children Colorado hig quality SNPs'] = params.chco_highqual_snps
    if (params.cadd_WG_SNVs_tbi)            summary['cadd_WG_SNVs_tbi']          = params.cadd_WG_SNVs_tbi
    

    if (workflow.profile == 'awsbatch') {
        summary['AWS Region']        = params.awsregion
        summary['AWS Queue']         = params.awsqueue
    }
    summary['Config Profile'] = workflow.profile
    if (params.config_profile_description)  summary['Config Description']  = params.config_profile_description
    if (params.config_profile_contact)      summary['Config Contact']      = params.config_profile_contact
    if (params.config_profile_url)          summary['Config URL']          = params.config_profile_url
    if (params.email) {
        summary['E-mail Address']        = params.email
        summary['MultiQC maxsize']       = params.maxMultiqcEmailFileSize
    }
    // params.properties.each { log.info "${it.key} -> ${it.value}" }
    log.info params.dump()
    log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
    if (params.monochrome_logs) log.info "----------------------------------------------------"
    else log.info "\033[2m----------------------------------------------------\033[0m"
    // log.info("ByeBye")
    // println(params)
}
def sortBedIntervalsByDescendingDuration(bedIntervals){
    bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
            final fields = line.split('\t')
            if (fields.size() >= 5) duration += fields[4].toFloat()
            else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotides_per_second
            }
        }
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    // println("passed list: $list")
    // println("real list: $realList")
    return list.every{ checkParameterExistence(it, realList) }
}

// Check if params.item exists and return params.genomes[params.genome].item otherwise
def checkParamReturnFile(item) {
    // Handle deprecation
    if (params.genomeDict && item == "dict") return file(params.genomeDict)
    if (params.genomeFile && item == "fasta") return file(params.genomeFile)
    if (params.genomeIndex && item == "fasta_fai") return file(params.genomeIndex)

    params."${item}" = params.genomes[params.genome]."${item}"
    return file(params."${item}")
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'fastqc',
        'samtools',
        'insert_size_metrics',
        'multiqc',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'mapping',
        'annotate',
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'starfusion',
        'arriba',
        'ericscript',
        'squid',
        'suppa'
    ]
}



// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}



// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
