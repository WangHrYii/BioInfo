nextflow.enable.dsl=2

params.reads = 'demo-data'
params.outdir = 'results/nextflow'

process COUNT_READS {
    tag "$sample"
    cpus 1
    memory '256 MB'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path summary_script

    output:
    path "${sample}.tsv"

    script:
    """
    python "${summary_script}" --sample "${sample}" \
      --r1 "${reads[0]}" --r2 "${reads[1]}" --out "${sample}.tsv"
    """
}

workflow {
    pairs = Channel.fromFilePairs("${params.reads}/*_R{1,2}.fastq.gz", checkIfExists: true)
    program = Channel.value(file("${projectDir}/../fastq_summary.py", checkIfExists: true))
    COUNT_READS(pairs, program)
}
