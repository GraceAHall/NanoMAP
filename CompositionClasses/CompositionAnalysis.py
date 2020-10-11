

from CompositionClasses.PafProcessing import PafProcessor
from CompositionClasses.Characterisation import Characteriser
from AbundanceClasses.ProportionEstimator import ProportionEstimator


import subprocess
import os


class CompositionAnalyser:
    def __init__(self, characterisation, context):
        self.characterisation = characterisation
        self.context = context


    def analyse(self):
        #self.set_genome_filenames()
        #self.prepare_database()
        #self.run_alignment()
        #self.process_paf()
        self.estimate_abundances()
        #self.create_final_characterisation()
        #self.write_characterisation()
        return self.identified_strains


    def estimate_abundances(self):
        pe = ProportionEstimator(self.characterisation, self.context)
        self.identified_strains = pe.estimate()


    def create_final_characterisation(self):
        ch = Characteriser(self.alignments, self.context)
        self.characterisation = ch.characterise()


    def write_characterisation(self):
        ct = self.context
        filename = ct.project_path + 'final_characterisation.tsv'
        self.characterisation.sort(key=lambda x: x.naive_abundance, reverse=True)
        print('\n\n{:>60}{:>60}{:>10}'.format('name', 'filename', 'sample abund.'))

        with open(filename, 'w') as fp:
            for s in self.characterisation:
                s.accessions.sort()
                accessions = '|'.join(s.accessions)

                print('{:>60}{:>60}{:>10.4f}'.format(s.name[-55:], s.filename, s.naive_abundance))
                
                for item in [s.name, s.filename, accessions, s.naive_abundance]:
                    fp.write(str(item) + '\t')
                fp.write('\n')



    '''
    def set_genome_filenames(self):
        folder = self.context.project_path + '/characterisations/'
        characterisation_filenames = os.listdir(folder)
        characterisation_filenames = [folder + x for x in characterisation_filenames]

        filenames = set()
        for filename in characterisation_filenames:
            print(filename)
            with open(filename, 'r') as fp:
                lines = fp.readlines()
                lines = [line.split('\t') for line in lines]

                try:
                    lines = [x[1] for x in lines if float(x[3]) > 0]
                except ValueError:
                    print(lines)

                for genomefile in lines:
                    filenames.add(genomefile)
    
        print(filenames)   
        self.genomes = list(filenames)


    def prepare_database(self):
        ct = self.context
        self.database_file = ct.project_path + '/identified_genomes.fa'

        if os.path.exists(self.database_file):
            os.remove(self.database_file)

        for genome in self.genomes:
            subprocess.call(f'cat {ct.database_path}{genome} >> {self.database_file}', shell=True)


    def estimate_abundances(self):
        pe = ProportionEstimator(self.characterisation, self.context)
        self.identified_strains_abundances = pe.estimate()


    def run_alignment(self):
        ct = self.context
        subprocess.call(f'minimap2 -c -x {ct.read_technology} -t {ct.threads} -I 1000G {self.database_file} {ct.fastq_path} > {ct.project_path}/composition.paf', shell=True)


    def process_paf(self):
        ct = self.context
        paf_file = ct.project_path + '/composition.paf'
        database_path = ct.database_path
        pp = PafProcessor(paf_file, database_path)
        alignments = pp.process()
        self.alignments = alignments


    

    '''