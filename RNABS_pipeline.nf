// Define input parameters

params.reads = "$baseDir/raw_data/test/*_{1,2}.fq"
params.GTF = "$baseDir/../Annotation/Mus_musculus.GRCm38.96.gtf"
params.GNM = "$baseDir/../Annotation/mm10.for.RNABS.fa"
params.multiqc = "$baseDir/multiqc_results"
params.outdir = "$baseDir/results"
params.rawdata = "$baseDir/raw_data"
params.working_data = "$baseDir/working_data"


log.info """\
         RNA Bisulfite Sequencing Pipeline    
         ===================================
         transcriptome: ${params.transcriptome_file}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         working_data : ${params.working_data}
         """
         .stripIndent()


// Check that paired reads exist, create read_id tuple: (pair_id, first file, second file)
Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .ifEmpty{ error "No matching reads"}
    .set { read_pairs_ch } 

process setup {

    """
    mkdir -p $params.outdir
    mkdir -p $params.working_data
    """
}

// Initial FastQC check
process fastqc {
    tag "FASTQC on $pair_id"
    publishDir "${params.rawdata}/${pair_id}", mode: 'copy'

    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch; read_pairs2_ch }
    


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
        .into { read_pairs_ch; read_pairs2_ch }


    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), file("*val*.fq") into clean_reads
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

    publishDir "${params.working_data}/${pair_id}",  mode: 'copy'

    input: 
    set val(pair_id), file(cleanReads) from clean_reads
    
    output:
    //stdout ch
    tuple pair_id, file("*.bam") into (raw_bam_Count, raw_bam_Call)
    set val(pair_id), file("*.bai") into raw_bai


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
    publishDir "${params.working_data}/${pair_id}",  mode: 'copy'

    input:
    set val(pair_id)

    output:
    set val(pair_id), file(cutoffFiles) into cutoff_bam

    """
    module reset
    module load Pysam
    python ${PWD}/bin/03b_readWrite_C_cutoff.py ${pair_id} ${PWD}
    """


}

// Call C sites
process callSites {
    tag "Counting Cs per read in $pair_id"
    scratch true
    cpus 40
    publishDir "${params.results}/${pair_id}",  mode: 'copy'

    input:
    set val(pair_id), file(mappedFile) from raw_bam_Call
    file(cutoffFiles) from cutoff_bam

    output:
    set val(pair_id), file(rawCountMatrix) into rawCounts
    
    script:
    """
    # Call raw map file
    module load meRanTK
    meRanCall -p 40 -bam ${mappedFile} -f ${params.GNM} -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10

    """

}



