// Define input parameters

params.reads = "$baseDir/raw_data/*/*_{1,2}.fq"
params.transcriptome_file = "$baseDir/../Annotation/mm10.fa"
params.multiqc = "$baseDir/multiqc_results"
params.outdir = "$baseDir/results"
params.rawdata = "$baseDir/raw_data"
params.working_data = "$baseDir/working_data"


log.info """\
         RNA Bisulfite Sequencing Pipeline    
         ===================================
         transcriptome: ${params.transcriptome_file}
         reads        : ${params.reads}
         raw_data     : ${params.rawdata}
         outdir       : ${params.outdir}
         working_data : ${params.working_data}
         """
         .stripIndent()


// Check that paired reads exist
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
    mkdir -p fastqc_${pair_id}_logs
    fastqc -f fastq -q ${reads} -o fastqc_${pair_id}_logs
    """  
}

// Clean raw reads
process cleanReads {
    tag "Cleaning $pair_id"
    scratch true
    
    publishDir "${params.working_data}/${pair_id}", pattern: '.fq.gz',  mode: 'copy'
    publishDir "${params.rawdata}/${pair_id}", pattern: '.json', mode: 'copy'
    publishDir "${params.rawdata}/${pair_id}", pattern: '.html', mode: 'copy'

    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch; read_pairs2_ch }


    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    file "*"

    script:
    """
    fastp -w ${task.cpus} -q 25 -f 6 -t 6 -l 50 --trim_poly_x --poly_x_min_len 10 \
     -i ${reads[0]} -I ${reads[1]} -o ${pair_id}_1.fq.gz -O ${pair_id}_2.fq.gz \
     --failed_out ${pair_id}_failed.fq.gz -j ${pair_id}.json -h ${pair_id}.html \
     --detect_adapter_for_pe --overlap_diff_percent_limit 25
    """  
}


