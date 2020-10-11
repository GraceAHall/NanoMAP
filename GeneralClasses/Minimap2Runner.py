

import os
import subprocess



class Minimap2Runner:
    def __init__(self, context):
        self.context = context

    
    def initial_alignment(self):
        print('running minimap2 on full database')
        ct = self.context

        paf_output = ct.project_path + '/runtimefiles/pafs/full_alignment.paf'
        preset = ct.read_technology
        threads = ct.threads
        index = ct.database_path + 'database.mmi'
        reads = ct.fastq_path
        
        with open(paf_output, 'w') as outfile:
            subprocess.run(['minimap2', '-x', preset, '-t', threads, index, reads], stdout=outfile)




# to delete
    '''
    def run_clades(self):
        ct = self.context
        group_names = os.listdir(ct.project_path + '/runtimefiles/group_databases/')
        group_names = [x for x in group_names if x.endswith('.fa')]
        group_names = [x.rsplit('.', 1)[0] for x in group_names]
        for group in group_names:
            database = ct.project_path + '/runtimefiles/group_databases/' + group + '.fa' 
            fastq = ct.project_path + '/runtimefiles/group_fastqs/' + group + '.fq' 
            filename = ct.project_path + '/runtimefiles/pafs/' + group + '.paf'
            with open(filename, 'w') as outfile:
                subprocess.run(['minimap2', '-c', '-x', ct.read_technology, '-t', ct.threads, '-K', '100M', '-I', '1000G', database, fastq], stdout=outfile)
    '''











# TO DELETE
    '''
    def run_pairwise(self, reference_genome_filenames):
        ct = self.context
        
        paf_file_count = 0
        reference_genome_pair_filename = ct.project_path + '/resolve_databases/reference_pair.fa'
        for i in range(len(reference_genome_filenames) - 1):
            for j in range(i + 1, len(reference_genome_filenames)):
                paf_file_count += 1

                try: 
                    os.remove(reference_genome_pair_filename)
                except FileNotFoundError:
                    pass

                genome1 = ct.database_path + reference_genome_filenames[i]
                genome2 = ct.database_path + reference_genome_filenames[j]
                paf_filename = ct.project_path + '/resolve_pafs/strain_group_' + ct.group_id + '_pair_' + str(paf_file_count) + '.paf'

                subprocess.call(f"cat {genome1} >> {reference_genome_pair_filename}", shell=True)
                subprocess.call(f"cat {genome2} >> {reference_genome_pair_filename}", shell=True)
                with open(paf_filename, "w") as outfile:
                    subprocess.run(['minimap2', '-c', '-x', ct.read_technology, '-t', ct.threads, '-K', '100M', '-I', '1000G', reference_genome_pair_filename, ct.fastq_file], stdout=outfile)
    '''

