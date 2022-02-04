// Define input parameters

params.reads = "$baseDir/raw_data/*_{1,2}.fq"
params.transcriptome_file = "$baseDir/../Annotation/mm10.fa"
params.multiqc = "$baseDir/multiqc_results"
params.outdir = "$baseDir/results"


log.info """\
         R N A S E Q - N F   P I P E L I N E    
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


