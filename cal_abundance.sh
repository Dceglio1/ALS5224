#!/usr/bin/sh

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=20G
#SBATCH -t 8:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3

eval "$(conda shell.bash hook)"

conda activate mypy3

#REF_arg='/projects/ciwars/ChezLiz/CARD4.0.1.dmnd'
REF_arg='/projects/ciwars/ChezLiz/protein_fasta_protein_homolog_model.len.txt' ##comment either this out or the deeparg length file depending on what you used. 
#REFA='/projects/ciwars/databases/bacmet_len.txt'

#REF_16s='/projects/ciwars/databases/gg_13_5.len'
REF_rpob='/projects/ciwars/databases/RpoB.ref.len'

#samples=$(ls $src/projects/ciwars/avdarling/mergedreadsfinal/fastp/S2* | awk '{split($_,x,"_fastp.fastq.gz"); print x[1]}') #change path to where your files are located and awk line depending on how you named your files
#samples=$(ls /projects/ciwars/thomasbyrne/*_clean_merged.fastq | awk '{split($_,x,"_clean_merged.fastq"); print x[1]}') #for my subset
cd /projects/ciwars/ChezLiz/vsearch_output
samples=$(ls R*merged.fasta.gz | awk -F/ '{gsub(/_clean_merged.fasta.gz/, "", $NF); print $NF}' | sort | uniq)

cd /projects/ciwars/ChezLiz/diamond_output
for sample in $samples; do
    #sample=$(basename "$sample")
    #echo $sample
    printf "%s \n" ${sample}
    file=$(basename -- ${sample})
    printf "%s \n" ${file}
    python rp_abun.py -a ${sample}_arg_short.csv -r ${sample}_rpob.csv -la $REF_arg -lr $REF_rpob --rpob_identity 40 --db card -o /projects/ciwars/ChezLiz/diamond_output/Abundances/${sample}_abundance_new.csv
done
