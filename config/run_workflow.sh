#!/bin/bash
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=bambu_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=preskaa@mskcc.org
#SBATCH --output=slurm%j_snkmk.out

snakemake=/data1/shahs3/users/preskaa/shared_conda_envs/snakemake/bin/snakemake
pipeline_dir=$HOME/bambu-smk
${snakemake} \
    --snakefile ${pipeline_dir}/workflow/Snakefile \
    --profile ${pipeline_dir}/workflow/profiles/ \
    --configfile ${pipeline_dir}/config/config.yml \
    --use-singularity \
    --singularity-args "--bind /data1/shahs3:/data1/shahs3" \
    --dry-run