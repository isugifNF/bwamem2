# bwamem2

```
git clone https://github.com/isugifNF/bwamem2.git
cd bwamem2
nextflow run main.nf --help
```

or 

```
nextflow run isugifNF/bwamem2 -r main --help
```



```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [voluminous_hodgkin] - revision: 5d7eca7c21

   Usage:
   The typical command for running the pipeline is as follows:
   nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile singularity
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
    --threads               Threads per process [default:4 for local, 16 for slurm]
    --window                Window size passed to bedtools for gatk [default:100000]
    --queueSize             Maximum jobs to submit to slurm [default:20]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help
```
