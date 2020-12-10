

import os
import subprocess


  

class IndexBuilder:
    def __init__(self, context):
        self.context = context


    def build(self):
        self.set_fasta_extensions()
        self.print_user_message()
        self.concatenate_fastas()
        self.build_index()
        self.remove_concatenated_fasta()


    def set_fasta_extensions(self):
        fasta_extensions = ['fasta', 'fa', 'fna']
        filenames = os.listdir(self.context.database_path)
        filenames = [f for f in filenames if '.' in f]

        fasta_filenames = [f for f in filenames if f.rsplit('.', 1)[1] in fasta_extensions]
        self.num_fastas = len(fasta_filenames)

        folder_extensions = set([f.rsplit('.', 1)[1] for f in filenames])
        folder_fasta_extensions = [ext for ext in folder_extensions if ext in fasta_extensions]
        self.fasta_extensions = folder_fasta_extensions


    def print_user_message(self):
        print('\n2. creating metagenome')
        print(f'{self.num_fastas} files to concatenate')
        if len(self.fasta_extensions) > 100:
            print('this could take a while')


    def concatenate_fastas(self):
        ct = self.context

        if os.path.exists(ct.database_path + 'database.fastaaaa'):
            os.remove(ct.database_path + 'database.fastaaaa')
        
        for ext in self.fasta_extensions:
            subprocess.call(f"cat {ct.database_path}*.{ext} >> {ct.database_path}database.fastaaaa", shell=True)  

        print('metagenome created.')        


    def build_index(self):
        print('\n3. building metagenome index')
        ct = self.context

        if os.path.exists(ct.database_path + 'database.mmi'):
            os.remove(ct.database_path + 'database.mmi')
            
        subprocess.call(f'minimap2 --idx-no-seq -x {ct.read_technology_preset} -I {str(ct.max_memory)}G, -t {ct.threads} -d {ct.database_path}database.mmi {ct.database_path}database.fastaaaa', shell=True)
        print('index created.')
        

    def remove_concatenated_fasta(self):
        if os.path.exists(self.context.database_path + 'database.fastaaaa'):
            os.remove(self.context.database_path + 'database.fastaaaa')