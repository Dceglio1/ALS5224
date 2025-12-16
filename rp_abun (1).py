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
