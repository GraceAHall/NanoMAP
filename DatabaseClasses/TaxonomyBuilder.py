


import os
import json
import sys
from collections import defaultdict


class TaxonomyBuilder:
    def __init__(self, context):
        self.context = context


    def build(self):
        print('\n1. building taxonomy')
        self.set_fastas()

        ag_dict, af_dict = self.make_dicts()
        fa_dict = self.invert_dict(af_dict)

        self.save_item(ag_dict, 'accessions_genomes.json') 
        self.save_item(af_dict, 'accessions_filenames.json') 
        self.save_item(fa_dict, 'filenames_accessions.json') 
        print('taxonomy built')


    def set_fastas(self):
        extensions = ['fa', 'fna', 'fasta']
        filenames = os.listdir(self.context.database_path)
        filenames = [f for f in filenames if '.' in f]
        fastas = [f for f in filenames if f.rsplit('.', 1)[1] in extensions]
        self.fastas = fastas


    def make_dicts(self):
        fastas = self.fastas
        folder = self.context.database_path
        accessions_genomenames = {}
        accessions_filenames = {}
        
        processed_file_count = 0
        total_files = len(fastas)
        print(f'{total_files} files to process')

        for filename in fastas:
            processed_file_count += 1
            if processed_file_count % 100 == 0:
                print(f'processed {processed_file_count} of {total_files} files')
            with open(folder + filename, 'r') as fp:
                line = fp.readline()
                while line:
                    if line.startswith('>'):
                        try:
                            accession, genomename = line.split(' ', 1)
                            genomename = genomename.rstrip('\n')
                            accession = accession.lstrip('>').upper()

                        except ValueError:
                            line = line.rstrip('\n').lstrip('>')
                            accession = line.upper()
                            genomename = line
                            
                        accessions_genomenames[accession] = genomename
                        accessions_filenames[accession] = filename
                    line = fp.readline()
                            
        return accessions_genomenames, accessions_filenames

    
    def invert_dict(self, the_dict):
        inverted_dict = {}
        for k, v in the_dict.items():
            inverted_dict[v] = k
        return inverted_dict


    def save_item(self, item, filename):
        taxonomy_folder = self.context.database_path + 'taxonomy/'
        if not os.path.exists(taxonomy_folder):
            os.mkdir(taxonomy_folder)

        with open(taxonomy_folder + filename, "w") as fp:
            json.dump(item, fp)





# Using NCBI taxdump - For later maybe!
    '''
    def build_ncbi_taxonomy(self):
        print('collecting taxonomy information')
        self.set_structure_filenames()
        self.associate_accessions_taxids()
        self.save_item(self.at_dict, 'accessions_taxids.json') 

 
    def set_structure_filenames(self):
        folder = self.context.taxonomy_path
        filenames = os.listdir(folder)
        filenames = [f for f in filenames if '_assembly_report.txt' in f]
        self.structure_filenames = filenames


    def associate_accessions_taxids(self):
        path = self.context.taxonomy_path
        accessions_taxids = {}
        for filename in self.structure_filenames:
            filepath = path + filename
            taxid, accessions = self.get_structure_report_info(filepath)
            for acc in accessions:
                accessions_taxids[acc] = taxid
        self.at_dict = accessions_taxids

    
    def get_structure_report_info(self, filepath):
        with open(filepath, 'r') as fp:
            lines = fp.readlines()

        taxid = self.get_taxid(lines)
        accessions = self.get_accessions(lines)
        return taxid, accessions


    @staticmethod
    def get_taxid(lines):
        taxid_line = [l for l in lines if "# Taxid:" in l][0]
        taxid = taxid_line.split(' ')[-1].rstrip('\n')
        return taxid


    @staticmethod
    def get_accessions(lines):
        for i, line in enumerate(lines):
            if '# Sequence-Name' in line:
                accession_line_start = i + 1
                break

        accession_lines = lines[accession_line_start:]

        accessions = []
        for line in accession_lines:
            line = line.split('\t')
            genbank = line[4]
            refseq = line[6]
            accessions.append(genbank)
            accessions.append(refseq)
        
        return accessions
    '''          
    