

 
from StrainGroupingClasses.PafProcessing import PafProcessor
from StrainGroupingClasses.Similarity import SimilarityCalculator
from StrainGroupingClasses.Grouping import Grouper
from StrainGroupingClasses.Abundance import GroupAbundanceCalculator
from StrainGroupingClasses.ReadBinning import ReadBinner
from StrainGroupingClasses.Pruning import GroupPruner
from StrainGroupingClasses.ReferenceMaker import ReferenceMaker
from StrainGroupingClasses.FastqMaker import FastqMaker
import sys
import json
from modules.io import load_json
import time



class StrainGrouper:
    def __init__(self, context):
        self.context = context


    def group(self):
        self.process_paf()
        self.calculate_strain_similarity()
        self.group_strains()
        self.calculate_group_abundances()
        self.bin_reads()
        self.prune_groups()
        #self.write_group_characterisation()
        self.create_reference_databases()
        self.create_fastqs()
        return self.strain_groups
        

    def process_paf(self):
        pp = PafProcessor(self.context.paf, self.context.database_path)
        self.read_classifications, self.strain_basecounts = pp.process()


    def calculate_strain_similarity(self):
        sc = SimilarityCalculator(self.read_classifications)
        self.read_counts, self.co_occurance = sc.pairwise_occurance()
        print('done')


    def group_strains(self):
        sg = Grouper(self.read_counts, self.co_occurance)
        self.strain_groups, self.group_membership = sg.group()


    def calculate_group_abundances(self):
        gac = GroupAbundanceCalculator(self.strain_groups, self.group_membership, self.strain_basecounts)
        self.strain_groups = gac.calculate()


    def bin_reads(self):
        rb = ReadBinner(self.strain_groups, self.group_membership, self.read_classifications)
        self.strain_groups = rb.bin()


    def prune_groups(self):
        gp = GroupPruner(self.strain_groups)
        self.strain_groups = gp.prune()
        self.print_groups()
        

    def write_group_characterisation(self):
        filename = self.context.project_path + '/runtimefiles/characterisations/group_level_characterisation.tsv'
        with open(filename, 'w') as fp:
            for group in self.strain_groups:
                fp.write(f'{group.id}\t{group.abundance}\n')


    def create_reference_databases(self):
        rm = ReferenceMaker(self.strain_groups, self.context)
        rm.make()


    def create_fastqs(self):
        fm = FastqMaker(self.strain_groups, self.context)
        fm.make()


    def print_groups(self):
        fa_dict = load_json(f'{self.context.database_path}taxonomy/filenames_accessions.json')
        ag_dict = load_json(f'{self.context.database_path}taxonomy/accessions_genomes.json')
        for group in self.strain_groups:
            print(f'\n\ngroup {group.id}: {group.abundance}')
            for identifier in list(group.identifiers)[:10]:
                accession = fa_dict[identifier]
                genome = ag_dict[accession.upper()]
                print(genome)


