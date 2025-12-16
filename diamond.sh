#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=64G
#SBATCH -t 144:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3
source activate diamond_and_vsearch

cd /projects/ciwars/ChezLiz/vsearch_output

ARGDB='/projects/ciwars/databases/CARD/CARD4.0.1.dmnd'
RPOB_DB='/projects/ciwars/databases/RpoB.dmnd'

sample=$(ls R*merged.fasta.gz | awk -F/ '{gsub(/_clean_merged.fasta.gz/, "", $NF); print $NF}' | sort | uniq)

for samples in $sample; do
	echo ${samples}
	diamond blastx -e 1e-10 --id 80 -k 1 --threads 64 -d $ARGDB -q ${samples}_clean_merged.fasta.gz \
	 -o /projects/ciwars/ChezLiz/diamond_output/${samples}_arg_short.csv --outfmt 6
	diamond blastx -e 1e-10 --id 40 -k 1 --threads 50 -d $RPOB_DB -q ${samples}_clean_merged.fasta.gz \
	 -o /projects/ciwars/ChezLiz/diamond_output/${samples}_rpob.csv --outfmt 6
done







