#!/bin/bash
#SBATCH -J NeoantigenML
#SBATCH -n 200
#SBATCH -p xhacnormalb
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --comment=test
#SBATCH -o %j
#SBATCH -e %j

source /public/home/yujijun01/miniconda3/bin/activate /public/home/yujijun01/miniconda3/envs/r41
Rscript /public/home/yujijun01/01_NeoantigenML/main_running_v3/06_mutate_kmer_single_combind_perfect_result/04_Model_running.R
