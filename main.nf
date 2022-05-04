#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMsg() {
  log.info """
   Usage:
   The typical command for running the pipeline is as follows:
   nextflow run main.nf --reference GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile singularity
   Mandatory arguments:
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm, singularity, docker [default:local]
    --singularity_img       Singularity image if [-profile singularity] is set [default:'shub://aseetharam/gatk:latest']
    --docker_img            Docker image if [-profile docker] is set [default:'j23414/gatk4']
    --bwamem2_app           Link to bwamem2 executable [default: 'bwa-mem2']
    --samtools_app          Link to samtools executable [default: 'samtools']
   Optional other arguments:
    --splitby               Number of sequence in each split of fastq.gz
    --threads               Threads per process [default:4 for local, 16 for slurm]
    --window                Window size passed to bedtools for gatk [default:100000]
    --queueSize             Maximum jobs to submit to slurm [default:20]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help
"""
}

if(params.help){
  helpMsg()
  exit 0
}

process bwamem2_index {
  tag "${genome_fasta.simpleName}"
  label 'bwamem'
  publishDir "${params.outdir}"

  input:
  path(genome_fasta)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("${genome_fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $bwamem2_app index $genome_fasta
  """
}

process bwamem2_mem {
  tag "$readname"
  label 'bwamem'
  publishDir "${params.outdir}"

  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(readpairs)

  output: // reads_mapped_2_genome.bam
  path("${readpairs.getAt(0).simpleName}_mapped.bam")

  script:
  """
  #! /usr/bin/env bash
  PROC1=\$(((`nproc`-1) * 3/4 + 1))
  PROC2=\$(((`nproc`-1) * 1/4 + 1))
  mkdir tmp
  ${bwamem2_app} mem -t \${PROC1} ${genome_fasta} ${readpairs} |\
     ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${readpairs.getAt(0).simpleName}_mapped.bam
  """
}
// samtools view --threads 1 -bS -

workflow {
    /* Input channels */
    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
    if(params.splitby) {
      reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true, flat:true)
        | splitFastq(by:"$params.splitby", compress:true, file:true, pe:true)
        | map { n -> [n.getAt(0), [n.getAt(1), n.getAt(2)] ] }
    } else {
      reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
    }
    
    /* main method */
    genome_ch | bwamem2_index | combine(reads_ch) | bwamem2_mem
}
