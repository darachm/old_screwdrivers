#!/usr/bin/env nextflow

/////
///// Pipeline for tabulating barseq data to counts, using barnone
/////

// To run you need the nextflow program. You can use the one at
//       /home/dhm267/nextflow
// or you can get your own by running the super-sketch command:
//       curl -s https://get.nextflow.io | bash
// Then you run something like:
//       ~/nextflow run main.nf
//   depending on where those two files are.

/////
///// Here's the easily tweakable parameters:
/////

 // This is for the report at the end. It's not the SLURM reports.
email_address = 'dchmiller@gmail.com'
 // The input gzipped fastq. 
gzipped_fastq = Channel.fromPath("/scratch/cgsb/gencore/out/Gresham/2017-04-14_AVA3A/1/000000000-AVA3A_l01n01.3310000008a5cf.fastq.gz")
 // The sample index file
sample_index = Channel.fromPath("sampleBarcodesRobinson2014.txt")
 // The strain index file
strain_index = Channel.fromPath("nislow_revised.txt")
 // How many different mismatches to try
mismatches_to_try = [0,1,2,3]

 // This is making a directory for the outputs and reports
file("./out").mkdirs()
file("./reports").mkdirs()

// If you want to edit more, check out the nextflow documentation.
// I use the single quoted shell blocks because it makes it clear
// what's a !{nextflow_variable} and a ${shell_variable}. Most 
// documentation doesn't use that yet.

/////
///// The actual processes
/////

 // This process gunzips everything into plain fastq
process gunzip {
  input: file fgz from gzipped_fastq
  output: file "fastq" into fastq
  shell: 
    '''
    zcat !{fgz} > "fastq"
    '''
}

 // This one actually launches out to slurm the barnone runs, one for
 // each mismatch parameter
process barnone {
  publishDir "out", mode: 'copy'
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
      -p 500000 \
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

