// Bare-bones pipeline for the analysis of RNA BS-seq data, based on the pipeline used in Johnson et al., 2022 NAR
// Setup: raw_data directory containing subdirectories with paired end reads
// Steps:
//  1) Adapter trimming, read cleaning, and mbias trimming of raw reads.
//  2) Mapping to C2T converted mm10 genome using meRanGh.
//  3) Performing a C-cutoff of 3 to mapped reads.
//  4) Methylation calling of signal (C-cutoff BAM) and noise (full BAM) using meRanGh


// Implementation wish list:
//  0) Option to use user-defined sample list as input (sample name, group name)
//  1) Option to use Gini-coefficient to determine C-cutoff
//  2) Signal/Noise filtering of called sites
//  3) Aggregation of final sites into Count Matrix
//  4) Addition of R analysis code to visualize results

// Define input parameters
// Currently testing: Option to use Gini coeff

params.reads = "$baseDir/raw_data/*/*_{1,2}.fq"
params.GTF = "$baseDir/../Annotation/Mus_musculus.GRCm38.96.gtf"
params.GNM = "$baseDir/../Annotation/mm10.for.RNABS.fa"
params.multiqc = "$baseDir/multiqc_results"
params.results = "$baseDir/results"
params.rawdata = "$baseDir/raw_data"
params.working_data = "$baseDir/working_data"
params.cutoff = "3"


log.info """\
         RNA Bisulfite Sequencing Pipeline    
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.results}
         working_data : ${params.working_data}
         C-cutoff     : ${params.cutoff}
         """
         .stripIndent()


// Check that paired reads exist, create read_id tuple: (pair_id, first file, second file)
Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .ifEmpty{ error "No matching reads"}
    .set { read_pairs_ch } 

process setup {

    """
    mkdir -p $params.results
    mkdir -p $params.working_data
    """
}

// Initial FastQC check
process fastqc {
    tag "FASTQC on $pair_id"
    publishDir "${params.rawdata}/${pair_id}", mode: 'copy'

    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch }
    
    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    path "fastqc_${pair_id}_logs" into fastqc_ch

    script:
    """
    module load FastQC
    mkdir -p fastqc_${pair_id}_logs
    fastqc -f fastq -q ${reads} -o fastqc_${pair_id}_logs
    """  
}

// Clean raw reads
process cleanReads {
    tag "Cleaning $pair_id"
    scratch true
    
    publishDir "${params.working_data}/${pair_id}",  mode: 'copy'  //pattern: '.fq.gz',
    //publishDir "${params.rawdata}/${pair_id}", pattern: '.json', mode: 'copy'
    //publishDir "${params.rawdata}/${pair_id}", pattern: '.html', mode: 'copy'

    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch }


    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), path("*val*.fq") into clean_reads
    //tuple val(pair_id), file "*.{html, json}" into read_stats 

    script:
    """
    module load fastp
    fastp -w ${task.cpus} -q 25 -f 6 -t 6 -l 50 --trim_poly_x --poly_x_min_len 10 \
     -i ${reads[0]} -I ${reads[1]} -o ${pair_id}_val_1.fq -O ${pair_id}_val_2.fq \
     --failed_out ${pair_id}_failed.fq.gz -j ${pair_id}.json -h ${pair_id}.html \
     --detect_adapter_for_pe --overlap_diff_percent_limit 25
    """  
}


// map reads to C2T converted genome
process mapReads {
    tag "Mapping reads from $pair_id"
    cpus 20
    scratch true

    publishDir "${params.working_data}/${pair_id}",  mode: 'copy', overwrite: true

    input: 
    set val(pair_id), file(cleanReads) from clean_reads
    
    output:
    //stdout ch
    val(pair_id) into pair_id_name
    set val(pair_id), path("*.bam") into raw_bam
    path("*.bai") into raw_bai


    script:
    """ 
    module load meRanTK
    meRanGh align -un -f ${cleanReads[0]} -r ${cleanReads[0]} -t 20 -fmo -ds -S ${pair_id}_genomeMap.sam \
     -ds -MM -id /projects/epigenomicslab/Annotation/mm10_meRanGh/BSgenomeIDX \
     -GTF /projects/epigenomicslab/Annotation/mm10.ensGene.for.RNABS.gtf
    """

}

// Count Cs per read in mapped files
process countCs {
    tag "Counting Cs per read in $pair_id"
    scratch true
    cpus 4
    publishDir "${params.working_data}/${pair_id}",  mode: 'copy', overwrite: true

    input:
    val(pair_id) from pair_id_name

    output:
    tuple val(pair_id), path("*_Ccutoff_PE.bam") into cutoff_bam
    val(pair_id) into pair_id_cutoff

    """
    module reset
    module load Pysam
    python ${PWD}/bin/03b_readWrite_C_cutoff.py ${pair_id} ${PWD}
    """


}

// Call C sites from raw mapped reads
process callSites {
    tag "Counting Cs per read in $pair_id"
    scratch true
    cpus 40
    publishDir "${params.results}/${pair_id}",  mode: 'copy'

    input:
    set val(pair_id), path(mappedFile) from raw_bam

    output:
    set val(pair_id), path("*_Call.txt") into rawCounts
    
        
    script:
    """
    # Call raw map file
    module load meRanTK
    meRanCall -p 40 \
    -bam ${mappedFile} \
    -o ${pair_id}_Call.txt \
    -f ${params.GNM} \
    -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10
    """

}


// Call C sites from cutoff file
process callCutoffSites {
    tag "Counting Cs per read in $pair_id"
    scratch true
    cpus 40
    publishDir "${params.results}/${pair_id}",  mode: 'copy'

    input:
    val(pair_id) from pair_id_cutoff

    output:
    set val(pair_id), path("*_cutoffCall.txt") into cutoffCounts
        
    script:
    """
    # Call raw map file
    module load meRanTK
    meRanCall -p 40 \
    -bam ${params.working_data}/${pair_id}/${pair_id}_${params.cutoff}_Ccutoff_PE.bam \
    -o ${pair_id}_${params.cutoff}_cutoffCall.txt \
    -f ${params.GNM} \
    -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10
    """

}


