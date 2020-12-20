

import subprocess
import os
from collections import defaultdict

from modules.io import load_genome
from modules.io import load_json
from modules.io import write_fastq
from modules.io import write_json
from AbundanceClasses.PafClassifier import PafClassifier
from AbundanceClasses.TerryEM import TerryEM
import math



class ProportionEstimator:
    def __init__(self, characterisation, group_abundance, context):
        self.present_strains = characterisation
        self.group_abundance = group_abundance
        self.context = context
        self.include_accessory_dna = False
        self.read_bases_mb = 10 * 1000000  # 20mb worth of reads per genome
        self.read_length = 3000


    def estimate(self):
        if len(self.present_strains) > 1:
            self.perform_alignment_based_estimation()
        elif len(self.present_strains) == 1:
            self.set_simple_abundance()
        else:
            print('No strains detected above abundance threshold in group. Continuing for other groups')
        return self.present_strains


    def perform_alignment_based_estimation(self):
        print('\nmultiple identified strains in group')
        print('performing composition estimation')
        self.perform_setup()
        self.create_database()
        self.load_genomes()
        self.create_fragments()

        #print('estimating transition probabilities')
        self.align_sample_reads()
        self.align_simulated_reads()

        self.create_symbols()
        self.create_accessions_filenames_mapping()
        empirical_pooled_summary, empirical_organism_summary = self.get_read_classifications(self.fragment_paf_file, empirical=True)
        observed_counts = self.get_read_classifications(self.sample_paf_file)

        path = self.context.project_path + '/runtimefiles/abundance_estimation/'
        write_json(empirical_pooled_summary, path + 'empirical_pooled_summary.json')
        write_json(empirical_organism_summary, path + 'empirical_organism_summary.json')
        write_json(observed_counts, path + 'observed_counts.json')

        empirical_organism_summary = load_json(path + 'empirical_organism_summary.json')
        observed_counts = load_json(path + 'observed_counts.json')

        print('EM rounds start')

        self.estimate_proportions(empirical_organism_summary, observed_counts)
        self.update_sample_proportions()


    def perform_setup(self):
        ct = self.context

        if not os.path.isdir(ct.project_path + '/runtimefiles/abundance_estimation/'):
            os.mkdir(ct.project_path + '/runtimefiles/abundance_estimation/')

        self.database_file = ct.project_path + '/runtimefiles/abundance_estimation/database.fa'
        self.sample_reads = ct.fastq_path
        self.sample_paf_file = ct.project_path + '/runtimefiles/abundance_estimation/sample_read_alignments.paf'
        self.fragment_reads = ct.project_path + '/runtimefiles/abundance_estimation/fragments.fq'
        self.fragment_paf_file = ct.project_path + '/runtimefiles/abundance_estimation/fragment_read_alignments.paf'
        

    def create_database(self):
        print('\npreparing reference database')
        ct = self.context
        
        if os.path.exists(self.database_file):
            os.remove(self.database_file)
        
        for strain in self.present_strains:
            subprocess.call('cat ' + ct.database_path + strain.filename + ' >> ' + self.database_file, shell=True)


    def load_genomes(self):
        genomes = []
        for strain in self.present_strains:
            filename = self.context.database_path + strain.filename
            new_genome = load_genome(filename)
            genomes.append(new_genome)
        self.genomes = genomes


    def create_fragments(self):
        filenames_reads = {}
        for genome in self.genomes:
            filename = genome.filename
            if self.include_accessory_dna:
                sequence = genome.get_all_dna()
            else:
                seqobject = genome.get_longest_sequence()
                sequence = seqobject.sequence
            filenames_reads[filename] = self.get_reads_stride(sequence)
            
        if filenames_reads is not None:
            write_fastq(filenames_reads, self.fragment_reads)
    

    def get_reads_stride(self, sequence):
        reads = []
        num_reads = math.floor(self.read_bases_mb / self.read_length)
        stride = math.ceil(len(sequence) / num_reads)
        sequence = sequence + sequence[:100000]
        for i in range(num_reads):
            extract_location = i * stride
            if extract_location > len(sequence):
                break
            reads.append(sequence[extract_location: extract_location + self.read_length])
        return reads
        

    def align_sample_reads(self):
        ct = self.context

        with open(self.sample_paf_file, 'w') as outfile:
            print('running minimap2 on database')
            subprocess.run(['minimap2', '-c', '-x', ct.read_technology, '-t', ct.threads, '-K', '50M', '-I', str(ct.max_memory) + 'G', self.database_file, self.sample_reads], stdout=outfile)


    def align_simulated_reads(self):
        ct = self.context

        with open(self.fragment_paf_file, 'w') as outfile:
            print('running minimap2 on database')
            subprocess.run(['minimap2', '-c', '-x', ct.read_technology, '-t', ct.threads, '-K', '50M', '-I', str(ct.max_memory) + 'G', self.database_file, self.fragment_reads], stdout=outfile)


    def create_symbols(self):
        filenames = []
        for genome in self.genomes:
            filenames.append(genome.filename)

        symbol_mappings = defaultdict(str)
        genome_symbols = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789!@#$%^&*()'
        symbol_pointer = 0

        for filename in sorted(filenames):
            symbol_mappings[filename] = genome_symbols[symbol_pointer]
            symbol_pointer += 1
        
        self.symbol_mappings = symbol_mappings

    
    def create_accessions_filenames_mapping(self):
        accessions_filenames = {}
        for genome in self.genomes:
            accessions = genome.get_accessions()
            for acccession in accessions: 
                accessions_filenames[acccession] = genome.filename
        self.accessions_filenames = accessions_filenames


    def get_read_classifications(self, paf_file, empirical=False):
        pc = PafClassifier(paf_file, self.accessions_filenames, self.symbol_mappings, self.context)
        pc.classify(empirical=empirical)
        #pc.print_classification_stats()
        if empirical:
            return pc.pooled_summary_dict, pc.organisms_summary_dict
        else:
            return pc.pooled_summary_dict


    def estimate_proportions(self, empirical_organism_summary, observed_counts):
        tem = TerryEM(empirical_organism_summary, observed_counts)
        self.sample_proportions, self.data_ll = tem.iterate()

    
    def update_sample_proportions(self):
        for terry_symbol, proportion in self.sample_proportions.items():
            #print(terry_symbol, proportion)
            for filename, mapping_symbol in self.symbol_mappings.items():               
                if terry_symbol == mapping_symbol:
                    for strain in self.present_strains:
                        if filename == strain.filename:
                            strain.sample_abundance = proportion * self.group_abundance
                            #print('{:>50}{:>15.2f}'.format(strain.name[:50], strain.sample_abundance))


    def set_simple_abundance(self):
        print('\nSingle strain detected in group. Setting sample abundance')
        self.present_strains[0].sample_abundance = self.group_abundance
        print()




