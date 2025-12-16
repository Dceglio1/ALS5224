# Using R to work with DIAMOND files

### Turning DIAMOND outputs into a workable file on ARC
DIAMOND outputs are a little tricky, as it does not compile the number of times a gene is seen in a sample. Therefore, we have to do it ourselves.

#### Code in ARC for DIAMOND

**diamond.sh**

Inputs:
- CARD4.0.1.dmnd: This is the DIAMOND file created from the fasta file from the reference database (CARD database, version 4.0.1)
- RpoB.dmnd: DIAMOND file showing all rpob genes. Created by someone else
- FASTQ files: Samples from after running fastq, bbduk, and vsearch

Note: --id only reports alignments above the listed percentage. 80% is chosen for the ARGs because the CARD database
is manually curated, and we are confident reads close to the reference reads are ARGs. 40% is chosen for the rpob
genes as those are much more common, and can have differences compared to the database. 

<details>
  <summary>diamond.sh</summary>
  
``` 
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

ARGDB='/projects/ciwars/ChezLiz/CARD4.0.1.dmnd'
RPOB_DB='/projects/ciwars/databases/RpoB.dmnd'

sample=$(ls R*merged.fasta.gz | awk -F/ '{gsub(/_clean_merged.fasta.gz/, "", $NF); print $NF}' | sort | uniq)

for samples in $sample; do
	echo ${samples}
	diamond blastx -e 1e-10 --id 80 -k 1 --threads 50 -d $ARGDB -q ${samples}_clean_merged.fasta.gz -o /projects/ciwars/ChezLiz/diamond_output/${samples}_arg_short.csv --outfmt 6
	diamond blastx -e 1e-10 --id 40 -k 1 --threads 50 -d $RPOB_DB -q ${samples}_clean_merged.fasta.gz -o /projects/ciwars/ChezLiz/diamond_output/${samples}_rpob.csv --outfmt 6
done

```
</details>

Outputs:
- 1 file for ARG alignments for each sample
- 1 file for rpob alignments for each sample
  
#### Making DIAMOND Outputs Usable
Next, put rp_abun.py in the same folder as your DIAMOND outputs.

This code was written by Monjura Afrin Rumi a few years ago. 
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


Then, run cal_abundances.sh. You'll first need to create an environment named "mypy3" and install python (version 3.1.2) and pandas (version 2.3.3) in the ARC. To do so, use the following lines of code.

```
conda create -n mypy3 python=3.12 pip 
source activate mypy3
conda install ipykernel
pip install plotly kaleido
conda install pandas
```
**cal_abundances.sh**

Run cal_abundances.sh to calculate the abundance of each ARG and RPOB gene in every sample. This code came from Thomas Byrne. 

Inputs:
- protein_fasta_protein_homolog_model.len.txt: This gives the length of each read in the referenee database. Made for CARD 4.0.1
- RpoB.ref.len: Length of each rpob gene in the reference database. File created by someone else (maybe Amanda Darling?)
- DIAMOND outputs: This includes the csv files for both ARGs and rpob genes for each sample

**Mistake to avoid: Make sure the --db is set to "card". The default is "deeparg", and without specifying the code will not run**

**Mistake to Avoid: Mak"e sure to include "eval "$(conda shell.bash hook)"", as this forces python to activate. I don't know why it doesn't without this...**
<details>
<summary> cal_abundances.sh</summary>

```
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

```
</details>

Outputs:
- 1 "abundance" file for each sample, showing rpob normalized ARG counts

**merging_my_normalized.sh**

Finally, run merging_my_noramlized.sh. This script merges the previous abundance_new.csv files for each sample into a single, final output: I2GDS_G6_AMR_Diamond.txt. This output will be one of the inputs for the R script.

Inputs:
- The "abundance" files for each sample from the previous step

<details>
<summary> merging_my_noramlized.sh</summary>

```
#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=20G
#SBATCH -t 8:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu
# Output file name

output_file="ChezLiz_merged_args_new.txt"

cd /projects/ciwars/ChezLiz/diamond_output/Abundances

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

Output:
- 1 merged file with all data, named "ChezLiz_merged_args_new.txt"

### Workin in R

**I2GDSFinal.R**

Below is the R code used. 

Packages: Tidyverse, version 2.0.0

**Mistake to Avoid: Make sure only TidyVerse is loaded. Other packages use the 'group_by' function, and will mess the code up if they're loaded**

Inputs:
- CARD4.0.1_aro_cat.csv: From the CARD database, contains metadata about each ARG. Received when downloading CARD
- ChezLiz_merged_args_new.txt: File from the previous steps
- ChezLizMetadata.csv: Metadata about the samples
- Antibiotic to Drug.csv: A csv I created from the deeparg database. Used to place ARGs into broader classes


<details>
<summary> I2GDSFinal.R</summary>

```
#Initializing tidyverse, setting working directory, inputting needed files
library(tidyverse)

setwd("C:/Users/Jvlan/Downloads")

#The drugs object comes from the deeparg database, where they assigned each specific antibiotic to a broader class
#The CARD_aro object is from the CARD database, and shows which gene confers resistance to each specific antibiotic
CARD_aro <- read.csv("CARD4.0.1_aro_cat.csv")
Args <- read.delim("ChezLiz_merged_args_new.txt")
metadata <- read.csv("ChezLizMetadata.csv")
drugs <- read.csv("Antibiotic to Drug.csv")

#To make the graph look better, we assign the genes in the CARD database to a broader drug class using the deeparg conversion
for (i in 1:nrow(CARD_aro)) {
  for (j in 1:nrow(drugs)) {
    if(grepl(drugs[j,1], CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- drugs[j,2]}
  }
}

#If a gene confers resistance to multiple drugs, we assign it "Multidrug"
for (i in 1:nrow(CARD_aro)) {
  if(grepl(";", CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- "Multidrug"}
}

#If it is still unassigned, it gets classified as other
CARD_aro$Drug[is.na(CARD_aro$Drug)] <- "Other"

#Merging the DIAMOND data and the CARD metadata about each gene
Args <- merge(Args, CARD_aro, by.x = "protein_accession", by.y = "Protein.Accession")

#Renaming the samples to only include the same name
Args$sample = substr(Args$sample, 1, nchar(Args$sample)-10)

#Merging the DIAMOND data and sample metadata
Args <- merge(Args, metadata, by.x = "sample", by.y = "sample_original")

#Just incase the rpob normalization messed up and left this as an NA
#Something went wrong in ARC if this is an NA!
Args$rpob_Normalization[Args$rpob_Normalization==""] <- 0

#Assigning data its appropriate type
Args$Drug.Class <- as.factor(Args$Drug.Class)
Args$Fraction <- as.factor(Args$Fraction)
Args$rpob_Normalization <- as.numeric(Args$rpob_Normalization)
Args$count <- as.numeric(Args$count)
Args$Drug <- as.factor(Args$Drug)

#Creating df where bubble plot will be generated from
#Fraction = Sampling location
#Only separating by fraction to see how different sampling locations are different
bubble <- Args %>%
  group_by(Drug, Fraction) %>%
  summarize(sum_rpob_normalized_count = mean(rpob_Normalization))

#The "CONTROL" sample type occurs every 5th row in bubble
controlRow <- 5

#Creating new column in bubble to put log differences
bubble$logDiff <- 0

#This loop goes through the df and calculates the log difference in 
#rpob normalized ARGs compared to the control for each drug class
for (i in 1:13) { #Number of drug classes
  for (k in 1:7) { #Number of sample types
    count = k+((i*7)-7)
    bubble[count, 4] <- log10(bubble[count, 3]/bubble[controlRow, 3])
  }
  controlRow <- controlRow +7
}

#Find absolute value for plotting
bubble$absLogDiff <- abs(bubble$logDiff)

#Differentiate between increased ARGs (logDiff > 0) and decreased (logDiff < 0)
bubble$Pos <- ifelse(bubble$logDiff >= 0, "Increase", "Decrease")

#Set order of Fraction for graph
bubble$Fraction <- factor(bubble$Fraction, levels = c("INF", "EFF", "BOIL", "30M", "100M", "Blank", "CONTROL"))

#Using ggplot to plot the graph
ggplot(subset(bubble, Fraction %in% c("100M", "30M", "BOIL", "EFF", "INF")), aes(x = Fraction, y = Drug, size = absLogDiff)) +
  geom_point(alpha = 0.8, aes(color = Pos)) +
  theme_bw() + 
  scale_size(range = c(1, 20)) +
  xlab("Sampling Location") +
  labs(size = "Log Difference to Control") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(color = "Change") 
```

</details>

Output:
- The image below

<img width="1020" height="723" alt="image" src="https://github.com/user-attachments/assets/05016fd1-02f9-4c3a-9754-0dacb709cdd2" />
Figure 1. Log change in average abundance of rpob gene normalzied ARGs in different sampling locations relative to the control location.  

The data comes from a wastewater treatment plant (WWTP) in Norfolk, VA, collected in 2021. INF and EFF are the influent and effluent of the WWTP respecively. After treatment, the effluent is sent into the Chesapeake Bay. BOIL represents water collected directly above the diffuser pipe (diffusing effluent into the Bay), while 30M and 100M are waters 30 meters and 100 meters setback from the diffuser. The Control was water collected at the Chesapeake Bay Bridge Tunnel, far away and (supposedly) free from influence from the WWTP.

Both INF and EFF have higher concentrations of all antibiotic resistance gene (ARG) classes compared to the Control. This is expected, as both WWTP influent and effluent are known to contain high concentrations of ARGs, relative to surface waters. Effluent waters had higher concentrations of sulfonamide ARGs compared to influent. Keenum et al. also found that sulfonamide and aminoglycoside ARGs increased from influent to effluent in chlorinated WWTPs. Effluent waters also showed increased concentrations of other drug and phenicol resistance encoding ARGs compared to influent, trends not observed in Keenum et al. All four 



Keenum, I., Calarco, J., Majeed, H., Hager-Soto, E. E., Bott, C., Garner, E., Harwood, V. J., & Pruden, A. (2024). To what extent do water reuse treatments reduce antibiotic resistance indicators? A comparison of two full-scale systems. Water Research, 254, 121425. https://doi.org/10.1016/j.watres.2024.121425

