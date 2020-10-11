

from collections import defaultdict


class Strain:
    def __init__(self, filename):
        self.filename = filename
        self.contigs = {}
        self.cec_ratio = 0
        self.name = ''
        self.basecount = 0
        self.naive_abundance = 0
        self.collinearity = 0
        self.mapq_alignments = []
        self.mapq_dict = defaultdict(int)
        self.group_id = None
        self.group_abundance = 0


    def set_basecount(self):
        self.basecount = sum([x.basecount for x in self.contigs.values()])


    def set_chromosomal_extrachromosomal_ratio(self):
        chromosome_basecount = sum([x.basecount for x in self.contigs.values() if x.is_chromosome])
        extra_chromosome_basecount = sum([x.basecount for x in self.contigs.values() if not x.is_chromosome])
        if extra_chromosome_basecount > 0:
            self.cec_ratio = chromosome_basecount / extra_chromosome_basecount
        else:
            self.cec_ratio = 999



class Contig:
    def __init__(self, accession, name):
        self.accession = accession
        self.name = name
        self.is_chromosome = True if 'plasmid' not in name and 'mitochondr' not in name else False
        self.basecount = 0




