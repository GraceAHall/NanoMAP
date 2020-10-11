


import subprocess
from GroupClasses.PafProcessor import PafProcessorLite
import os

class SubsetAligner:
    def __init__(self, filenames, group_reads, context):
        self.filenames = filenames
        self.group_reads = group_reads
        self.context = context


    def run(self):
        self.concatenate_genomes()
        self.align()
        self.characterise()
        return self.characterisation


    def concatenate_genomes(self):
        outputfile = 'runtimefiles/genome_subset.fasta'
        if os.path.exists(outputfile):
            os.remove(outputfile)
        for filename in self.filenames:
            subprocess.call('cat ' + self.context.database_folder + filename + ' >> ' + outputfile, shell=True)


    def align(self):
        outputfile = 'runtimefiles/genome_subset.paf'
        with open(outputfile, 'w') as outfile:
            subprocess.run(['minimap2', '-c', '-x', 'map-pb', '-t', self.context.threads, 'runtimefiles/genome_subset.fasta', self.group_reads], stdout=outfile)

#
    def characterise(self):
        ppl = PafProcessorLite('runtimefiles/genome_subset.paf')
        return ppl.process()






