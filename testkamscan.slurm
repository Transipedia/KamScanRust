#!/bin/bash
#SBATCH --job-name="kamscan"
#SBATCH -t 24:00:00 --mem 20G
#SBATCH --cpus-per-task=2


KAM="target/release/kamscan"
TSV="test/test.mat.tsv"
DESIGN="test/test.design"
OUT="test/test.out"
THREADS=2

time $KAM --tsv $TSV -d $DESIGN -o $OUT



