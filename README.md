# Working with data in R to create a stacked bar plot

### Turning DIAMOND outputs into a workable file
DIAMOND outputs are a little tricky, as it does not compile the number of times a gene is seen in a sample. Therefore, we have to do it ourselves.

First, run DIAMOND again, this time annotating against an rpob reference database. The necessary reference file (RpoB.dmnd) is located in code-files.

<details>
  <summary>diamondrpob.sh</summary>
  
``` 
#!/bin/bash
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
#SBATCH --cpus-per-task=16

module load Miniconda3

source activate diamond_and_vsearch

RPOB_DB='/projects/ciwars/databases/RpoB.dmnd' #replace with location of your downloaded RpoB.dmnd file

cd #insert path to directory where all vsearch output files are located

samples=$(ls *_clean_merged.fastq | awk '{gsub(/_clean_merged.fastq/,"",$0); print}' | sort | uniq)

for sample in $samples
do
	diamond blastx -e 1e-10 --id 40 -k 1 --threads 10 -d $RPOB_DB -q ${sample}_clean_merged.fastq -o ${sample}_rpob.csv --outfmt 6
done
```
</details>

Next, put rp_abun.py in the same folder as your DIAMOND outputs.
<details>
  <summary>rp_abun.py</summary>
  
``` 
import os
import argparse
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=Warning)

class annotated(object):
    def __init__(self, param):
        self.file = param['filename']
        self.length_file = param['length_file']
        self.identity = float(param['identity'])
        self.mlen = int(param['mlen'])
        self.evalue = float(param['evalue'])
 

    def data_process(self):
        data = pd.read_csv(self.file, sep = '\t', header = None) 
        data.sort_values([0,11], inplace = True, ascending = False)
        data.drop_duplicates(subset = 0, keep = 'first', inplace = True)
        data.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 
                'sEnd', 'evalue', 'bit']

        data = data[data.identity >= self.identity]
        data = data[data.alignLen >= self.mlen]
        data = data[data.evalue <= self.evalue]
        return data
    
    def length_process(self, database):
        len_data = pd.read_csv(self.length_file, sep = '\t', header = None) 
        len_data.columns = ['sub_id', 'gene_len']
        len_temp = len_data['sub_id'].str.split('|', expand=True)
        len_temp.columns = database.colum_names
        len_temp = len_temp.drop(columns = database.dropped_column)
        len_temp['gene_len'] = len_data.iloc[:, 1]
        len_temp['gene_len(bp)'] = (len_temp['gene_len']+1)*3
        return len_temp
    
    def combine_length(self, data, len_data, database):
        temp = data['sub_id'].str.split('|', expand=True)
        temp.columns = database.colum_names
        temp = temp.drop(columns = database.dropped_column)
        temp['length'] = data.iloc[:, 3]
        data = temp.groupby(database.keep_column)['length'].agg([('count','count'), ('length','sum')]).reset_index()
        data = pd.merge(data, len_data, how = "left", on = database.keep_column)
        return data
        

class Database(object):
    def __init__(self, name):
        self.colum_names = []
        self.dropped_column = []
        self.keep_column = []
        self.set_properties(name)
        
    def set_properties(self, name):
        if name == 'deeparg':
            self.colum_names = ['protein_accession', 'extra', 'gene', 'drug', 'gene_family']
            self.dropped_column = ['extra']
            self.keep_column = ['gene', 'drug', 'protein_accession', 'gene_family']
        elif name == 'card':
            self.colum_names = ['extra','protein_accession', 'aro_index', 'gene']
            self.dropped_column = ['extra']
            self.keep_column = ['gene', 'protein_accession', 'aro_index']
        elif name == 'rpob':
            self.colum_names = ['extra','accession_num', 'rpob_gene']
            self.dropped_column = ['extra']
            self.keep_column = ['accession_num', 'rpob_gene']
        elif name == 'mobileOG':
            self.colum_names = ["mobileOG Entry Name", "col1", "col2", "mge_class", "mge_subclass", "col5", "col6"]
            self.dropped_column = ["col1", "col2", "col5", "col6"]
            self.keep_column = ['mobileOG Entry Name', 'mge_class', 'mge_subclass']


def cal_abundance(arg_param, rpob_param, database, output_file = ""):
    arg_db = Database(database)
    rpob_db = Database('rpob')
    arg_obj = annotated(arg_param)
    rpob_obj = annotated(rpob_param)
    
    arg_data = arg_obj.data_process()
    arg_len_data = arg_obj.length_process(arg_db)
    arg_data = arg_obj.combine_length(arg_data, arg_len_data, arg_db)
    rpob_data = rpob_obj.data_process()
    rpob_len_data = rpob_obj.length_process(rpob_db)
    rpob_data = arg_obj.combine_length(rpob_data, rpob_len_data, rpob_db)
    arg_db.keep_column.append("count")
    output = arg_data[arg_db.keep_column]

    Nrpob = rpob_data.loc[rpob_data['length'] >= 0]['count'].sum()
    Lrpob = rpob_len_data['gene_len(bp)'].mean()
    

    output['rpob_Normalization'] = (arg_data['count']/arg_data['gene_len(bp)'])/(Nrpob/Lrpob)
    output['rpob_count'] = Nrpob 
    output['sample'] = os.path.splitext(os.path.basename(arg_param['filename']))[0]
    if not output_file:
        output_file = os.path.join(os.path.dirname(arg_param['filename']), os.path.splitext(os.path.basename(arg_param['filename']))[0] + "_abundance.txt")
    
    
    output.to_csv(output_file, sep = "\t", index = False)

    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--arg", type = str, required = True, help = "Path to arg annotation file (required)")
    parser.add_argument("-la", "--len_arg", type = str, required = True, help = "Path to arg database length file (required)") 
    parser.add_argument("-r", "--rpob", type = str, required = True, help = "Path to rpob annotation file (required)")
    parser.add_argument("-lr", "--len_rpob", type = str, help = "Path to rpob database length file")
    parser.add_argument("-o", "--out", type = str, help = "Path to output file")

    parser.add_argument("--db", type = str, default = "deeparg", help = "name of ARG database such as card/deeparg [default deeparg]")
    parser.add_argument("--arg_identity", type = float, default = 80, help = "minimum identity for alignments [default 80]",)
    parser.add_argument("--arg_mlen", type = float, default = 25, help = "diamond minimum length for considering a hit [default 25aa]",)
    parser.add_argument("--arg_evalue", type = float, default = 1e-10, help = "minimum e-value for alignments [default 1e-10]",)

    parser.add_argument("--rpob_identity", type = float, default = 40, help = "minimum identity for alignments [default 40]",)
    parser.add_argument("--rpob_mlen", type = float, default = 25, help = "diamond minimum length for considering a hit [default 25aa]",)
    parser.add_argument("--rpob_evalue", type = float, default = 1e-10, help = "minimum e-value for alignments [default 1e-10]",)

    args = parser.parse_args()
    db_name = args.db
    arg_parameters = dict(
        filename = args.arg,
        length_file = args.len_arg, 
        identity = args.arg_identity,
        mlen = args.arg_mlen,
        evalue = args.arg_evalue
    )

    rpob_parameters = dict(
        filename = args.rpob,
        length_file = args.len_rpob,
        identity = args.rpob_identity,
        mlen = args.rpob_mlen,
        evalue = args.rpob_evalue
    )
    
    
    cal_abundance(arg_param = arg_parameters, rpob_param = rpob_parameters, database = db_name, output_file = args.out)

if __name__ == '__main__':
    main()
```
</details>


Then, run cal_abundances.sh. You'll first need to create an environment named "mypy3" and install python and pandas. To do so, use the following lines of code.
```
conda create -n mypy3 python=3.12 pip 
source activate mypy3
conda install ipykernel
pip install plotly kaleido
conda install pandas
```
Run cal_abundances.sh to calculate the abundance of each ARG and RPOB gene in every sample.
<details>
<summary> cal_abundances.sh</summary>

```
#!/usr/bin/sh

#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu

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
cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output
samples=$(ls S*_arg_full.csv | awk -F/ '{gsub(/_arg_full.csv/, "", $NF); print $NF}' | sort | uniq)

for sample in $samples; do
    #sample=$(basename "$sample")
    #echo $sample
    printf "%s \n" ${sample}
    file=$(basename -- ${sample})
    printf "%s \n" ${file}
    python rp_abun.py -a ${sample}_arg_full.csv -r ${sample}_rpob.csv -la $REF_arg -lr $REF_rpob --rpob_identity 40 --db card -o /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output/${sample}_abundance_new.csv
done
```
</details>

Finally, run merging_my_noramlized.sh. This script merges the previous abundance_new.csv files for each sample into a single, final output: I2GDS_G6_AMR_Diamond.txt. This output will be one of the inputs for the R script.
<details>
<summary> merging_my_noramlized.sh</summary>

```
#!/bin/bash

#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --mem=1G
#SBATCH -t 1:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
# Output file name

output_file="I2GDS_G6_AMR_Diamond.txt"

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output

# Find and merge files

for file in *_abundance_new.csv; do
    # Skip the output file to avoid self-merging
    if [ "$file" != "$output_file" ]; then

        # Print filename as a comment in the merged file. nope i made it a column so that i can just treat it as a csv directly
        #echo "# $file" >> "$output_file"

        # Append contents to the merged file, skipping the header line
        #tail -n +2 "$file" >> "$output_file"

        #so i think i just make it 1 instead of 2 but we'll see ig
        tail -n +1 "$file" >> "$output_file"
    fi
done


echo "Merge complete. Merged file: $output_file"
```

</details>
