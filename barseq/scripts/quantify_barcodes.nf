#!/usr/bin/env nextflow

/////
///// Pipeline for tabulating barseq sequencing data to counts
/////

email_address = 'dchmiller@gmail.com'
gzipped_fastq = Channel.fromPath("/scratch/cgsb/gencore/out/Gresham/2017-04-14_AVA3A/1/000000000-AVA3A_l01n01.3310000008a5cf.fastq.gz")
sample_index = Channel.fromPath("sample_barcodes_robinson2014.txt")
strain_index = Channel.fromPath("mutant_barcodes_smith2009revision.txt")
mismatches_to_try = [0,1,2,3]

file("./tmp").mkdirs()
file("./reports").mkdirs()

process gunzip {
  input: file fgz from gzipped_fastq
  output: file "fastq" into fastq
  shell: 
    '''
    zcat !{fgz} > "fastq"
    '''
}

process barnone {
  publishDir "tmp", mode: 'copy'
  input:
    file fastq
    file sample_index
    file strain_index
    each MM from mismatches_to_try
  output:
    set file("barnone_output_${MM}.counts"), 
      file("barnone_output_${MM}.mismatch"), 
      file("barnone_output_${MM}.revised"),
      file("barnone_output_${MM}.stdout") into results
  shell: 
    '''
    module purge
    module load barnone/intel/20170501
    BarNone -f fastq \
      --multiplexfile !{sample_index} \
      --multiplexstart 1 --multiplexlength 5 \
      --tagstart 18 --taglength 3 --start 21 \
      --mismatches !{MM} \
      --mismatchfile barnone_output_!{MM}.mismatch \
      --revisedcatalog barnone_output_!{MM}.revised \
      -p 100000 \
      !{fastq} \
      barnone_output_!{MM}.counts \
      !{strain_index} > barnone_output_!{MM}.stdout
    '''
}

 // Special trigger for `onComplete`. I copied this from documentation.
 // Some predefined variables. It somehow mails it. Cool.
workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

  def subject = 'barnone run'
  def recipient = email_address

  ['mail', '-s', subject, recipient].execute() << """

  Pipeline execution summary
  ---------------------------
  Completed at: ${workflow.complete}
  Duration    : ${workflow.duration}
  Success     : ${workflow.success}
  workDir     : ${workflow.workDir}
  exit status : ${workflow.exitStatus}
  Error report: ${workflow.errorReport ?: '-'}
  """
}

