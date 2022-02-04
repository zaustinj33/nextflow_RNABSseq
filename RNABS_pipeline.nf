// Define input parameters

params.reads = "$baseDir/raw_data/*_{1,2}.fq"
params.transcriptome_file = "$baseDir/../Annotation/mm10.fa"
params.multiqc = "$baseDir/multiqc_results"
params.outdir = "$baseDir/results"


log.info """\
         RNA Bisulfite Sequencing Pipeline    
         ===================================
         transcriptome: ${params.transcriptome_file}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         working_data : ${params.workingdata}
         """
         .stripIndent()


// Check that paired reads exist
Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

process setup {

    """
    mkdir $params.outdir
    mkdir $params.workingdata
    """
}


process fastqc {
    tag "FASTQC on $pair_id"
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch; read_pairs2_ch }

    input:
    tuple pair_id, path(reads) from read_pairs_ch

    output:
    path "fastqc_${pair_id}_logs" into fastqc_ch

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
    """  
}

process cleanReads {
    tag "Cleaning $pair_id"
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .into { read_pairs_ch; read_pairs2_ch }
    publishDir "$params.workingdata"

    input:
    tuple pair_id, path(reads) from read_pairs_ch

    output:
    path "fastqc_${pair_id}_logs" into fastqc_ch

    script:
    """
    fastp -w 16 -q 25 -f 6 -t 6 -l 50 --trim_poly_x --poly_x_min_len 10 -i${reads[0]} -I ${reads[1]} --out1 ${pair_id}_1_qual.fq.gz --out2 ${pair_id}_2_qual.fq.gz
     --failed_out ${pair_id}_failed.fq.gz -j ${pair_id}.json -h ${pair_id}.html --detect_adapter_for_pe --overlap_diff_percent_limit 25
    """  
}


