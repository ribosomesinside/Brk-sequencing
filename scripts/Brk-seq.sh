#! /usr/bin/bash
#SBATCH --job-name=ubi_brk
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --partition=cpu
#SBATCH -o test_results/run.o
#SBATCH -e test_results/run.e
#SBATCH --mail-user=maria.mikhaylova@crick.ac.uk
#SBATCH --mail-type=BEGIN,FAIL,END

ml purge
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml CAMP_proxy

export WORK=/nemo/lab/vincentj/home/users/mikhaym/02.nf-seq/work
export NXF_HOME=/nemo/lab/vincentj/home/users/mikhaym/02.nf-seq/nextflow
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/vincentj/home/users/mikhaym/02.nf-seq/singularity

nextflow run nf-core/rnaseq -profile crick \
    --input /nemo/lab/prieto-godinoll/working/Jason/RNAseq/dros_samples/dros_test5/samplesheet/samplesheet.csv \
    --outdir /nemo/lab/prieto-godinoll/working/Jason/RNAseq/dros_samples/dros_test6/ \
    --max_memory 64GB \
    --aligner star_rsem \
    --fasta /nemo/lab/prieto-godinoll/working/Jason/RNAseq/dros_samples/genomes/Dmel/fasta/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa \
    --gtf /nemo/lab/prieto-godinoll/working/Jason/RNAseq/dros_samples/genomes/Dmel/gtf/Drosophila_melanogaster.BDGP6.32.106.gtf \
    --save_reference \
    --skip_biotype_qc \
    --email jason.somers@crick.ac.uk \
    -r 3.10 \
    -resume