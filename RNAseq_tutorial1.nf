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
         """
         .stripIndent()


// Check that paired reads exist
Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

process setup {

    """
    mkdir $params.outdir
    """
}

process index {
    input:
    path transcriptome from params.transcriptome_file
    
    output:
    path 'index' into index_ch
    
    script:
    println "task cpus: $task.cpus"
    """
    salmon index --threads 8 -t $transcriptome -i index
    """
}


process quantification {
    
    tag "Salmon on $pair_id"
    publishDir params.outdir, mode:'copy'

    input:
    path index from index_ch
    tuple pair_id, path(reads) from read_pairs_ch
 
    output:
    path pair_id into quant_ch
 
    script:
    """
    salmon quant --threads 8 --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
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

