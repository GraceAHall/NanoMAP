


from GroupClasses.PafProcessing import PafProcessor
from GroupClasses.Characterisation import Characteriser
from GroupClasses.Completion import CompletionJudge
from GroupClasses.Characterisation import Halver
from GroupClasses.StrainPicking import StrainPicker
from AbundanceClasses.ProportionEstimator import ProportionEstimator

import sys
import os
import subprocess



class GroupAnalyser:
    def __init__(self, group, context):
        self.group = group
        self.context = context
        self.include_plasmids_mitochondria = context.include_plasmids_mitochondria
        self.fastq = context.project_path + f'/runtimefiles/group_fastqs/strain_group_{group.id}.fq'
        self.paf = context.project_path + f'/runtimefiles/pafs/strain_group_{group.id}.paf'
        self.complete_database = context.project_path + f'/runtimefiles/group_databases/strain_group_{group.id}_initial.fa'
        self.narrow_database = context.project_path + f'/runtimefiles/group_databases/strain_group_{group.id}_narrow.fa'
        self.resolved = False


    def analyse_group(self):
        print('\nnarrowing strain candidates:\n')
        self.perform_alignment(initial=True)
        self.process_paf()
        self.create_characterisation()
        self.print_characterisation()
        self.judge_completion()
       
        while not self.resolved:
            self.halve_candidates()
            self.update_database()
            self.perform_alignment()
            self.process_paf()
            self.create_characterisation()
            self.print_characterisation()
            self.judge_completion()

        self.write_to_detailed_report()
        self.pick_strains()
        self.print_identified_strains()
        self.estimate_abundances()
        print(f'group {self.group.id} narrowing done.')

        return self.identified_strains



    def perform_alignment(self, initial=False):
        ct = self.context
        database = self.complete_database if initial else self.narrow_database

        with open(self.paf, 'w') as outfile:
            subprocess.run(['minimap2', '-c', '-p', '0.1', '--no-long-join', '-N', '10', '-x', ct.read_technology, '-I', str(ct.max_memory) + 'G', '-t', ct.threads, '-I', '1000G', database, self.fastq], stdout=outfile)
        

    def process_paf(self):
        pp = PafProcessor(self.paf, self.context.database_path, self.include_plasmids_mitochondria)
        self.alignments, self.collinearities = pp.process()


    def create_characterisation(self):
        ch = Characteriser(self.alignments, self.collinearities, self.group.id, self.group.abundance, self.context)
        self.characterisation = ch.format_characterisation()


    def judge_completion(self):
        cj = CompletionJudge(self.characterisation)
        self.resolved = cj.judge()


    def halve_candidates(self):
        hv = Halver(self.characterisation, self.context)
        self.characterisation = hv.halve()


    def update_database(self):
        filenames = [strain.filename for strain in self.characterisation]
        filenames = [self.context.database_path + f for f in filenames]
        
        try:
            os.remove(self.narrow_database)
        except FileNotFoundError:
            pass
        
        for f in filenames:
            subprocess.call(f'cat {f} >> {self.narrow_database}', shell=True)


    def write_to_detailed_report(self):
        output_filepath = self.context.project_name + '_detailed_report.tsv'
        self.characterisation.sort(key=lambda x: x.mapq_dict[60], reverse=True)
        for strain in self.characterisation:
            with open(output_filepath, 'a') as fp:
                fp.write(f'{strain.name}\t{strain.filename}\t{self.group.id}\t{round(strain.naive_abundance, 2)}\t{strain.mapq_dict[60]}\t{strain.mapq_dict[10]}\t{strain.mapq_dict[2]}\n')


    def pick_strains(self):
        sp = StrainPicker(self.characterisation)
        self.identified_strains = sp.pick()


    def estimate_abundances(self):
        pe = ProportionEstimator(self.identified_strains, self.group.abundance, self.context)
        self.identified_strains = pe.estimate()


    def print_identified_strains(self):
        print('\nidentified strains:')
        for strain in self.identified_strains:
            print(strain.name)
       

    def print_characterisation(self):
        self.characterisation.sort(key=lambda x: x.naive_abundance, reverse=True)
        print('\n\n{:>60}{:>60}{:>15}{:>10}{:>10}{:>10}{:>10}'.format('name', 'filename', 'naive abund.', 'MAPQ_60', 'MAPQ_10', 'MAPQ_2', 'MAPQ_1'))
        for s in self.characterisation:
            print('{:>60}{:>60}{:>15.2f}{:>10}{:>10}{:>10}{:>10}'.format(s.name[-55:], s.filename, s.naive_abundance, s.mapq_dict[60], s.mapq_dict[10], s.mapq_dict[2], s.mapq_dict[1]))


    def report_characterisation(self):
        filename = self.context.project_path + '/runtimefiles/characterisations/group_' + self.context.group_id + '.tsv'
        self.characterisation.sort(key=lambda x: x.sample_abundance, reverse=True)
        print('\n\n{:>60}{:>60}{:>15}{:>10}{:>10}{:>10}{:>10}'.format('name', 'filename', 'sample abund.', 'MAPQ_60', 'MAPQ_10', 'MAPQ_2', 'MAPQ_1'))

        with open(filename, 'w') as fp:
            for s in self.characterisation:
                print('{:>60}{:>60}{:>15.2f}{:>10}{:>10}{:>10}{:>10}'.format(s.name[-55:], s.filename, s.sample_abundance, s.mapq_dict[60], s.mapq_dict[10], s.mapq_dict[2], s.mapq_dict[1]))
                s.accessions.sort()
                accessions = '|'.join(s.accessions)
                for item in [s.name, s.filename, accessions, s.sample_abundance]:
                    fp.write(str(item) + '\t')
                fp.write('\n')


    






'''
    def perform_setup(context):
        if not os.path.exists(context.project_path + '/runtimefiles/'):
            os.mkdir(context.project_path + '/runtimefiles/')
'''