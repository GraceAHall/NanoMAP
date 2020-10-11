

from collections import defaultdict



class StrainGroup: 
    def __init__(self, group_id):
        self.group_id = group_id
        self.strains = []
        self.read_ids = set()
        self.filenames = set()
        #self.identifiers = []
        self.sample_abundance = 0
        self.is_sample_candidate = False
        self.mapqs = []


    def add(self, new_strain):
        self.sample_abundance += new_strain.group_abundance
        self.read_ids = self.read_ids.union(new_strain.read_ids)
        self.filenames.add(new_strain.filename)
        self.strains.append(new_strain)
    

    def get_cumulative_strain_abundance(self):
        cum_strain_abundance = 0
        for strain in self.strains:
            cum_strain_abundance += strain.group_abundance
        return cum_strain_abundance







class Strain:
    def __init__(self, filename):
        self.filename = filename
        self.is_sample_candidate = False
        self.mapq_confirmed = False
        self.mapq_60_count = 0
        self.mapq_10_count = 0
        self.group_abundance = 0
        self.sample_abundance = 0
        self.full_blocks = 0
        self.truncated_blocks = 0
        self.sequences = {} # starts as dict, gets reformatted. Maybe should have 2 classes actually. 
 
    
    def get_accessions(self):
        return [x.accession for x in self.sequences]


    def get_names(self):
        return [x.name for x in self.sequences]


    def get_primary_name(self):
        chromosomal_sequences = self.get_chromosome_sequences()
        chromosome_names = [x.name for x in chromosomal_sequences]
        chromosome_names.sort()
        primary_name = chromosome_names[0]
        primary_name = primary_name.rsplit('chromosome', 1)[0]
        primary_name = primary_name.rsplit(', complete', 1)[0]
        return primary_name


    def get_chromosome_read_ids(self):
        if not hasattr(self, 'chromosome_read_ids'):
            chromosomal_sequences = self.get_chromosome_sequences()
            read_ids = []
            for sequence in chromosomal_sequences:
                read_ids += sequence.get_read_ids()
            self.chromosome_read_ids = set(read_ids)
        return self.chromosome_read_ids
    

    def get_all_read_ids(self):
        if not hasattr(self, 'all_read_ids'):
            read_ids = []
            for sequence in self.sequences:
                read_ids += sequence.get_read_ids()
            self.read_ids = set(read_ids)
        return self.read_ids


    def get_total_basecount(self):
        if not hasattr(self, 'total_basecount'):
            self.total_basecount = sum([x.get_basecount() for x in self.sequences])
        return self.total_basecount

    
    def get_hma_counts(self):
        chromosomal_sequences = self.get_chromosome_sequences()
        self.mapq_60_count = 0
        self.mapq_10_count = 0

        for sequence in chromosomal_sequences:
            m60_count, m10_count = sequence.get_hma_count()
            self.mapq_60_count += m60_count
            self.mapq_10_count += m10_count

        return self.mapq_60_count, self.mapq_10_count


    def get_chromosomal_extrachromosomal_ratio(self):
        chromosomal_sequences = self.get_chromosome_sequences()
        extrachromosomal_sequences = self.get_extrachromosomal_sequences()

        chromosomal_basecount = 10
        extrachromosomal_basecount = 1

        for sequence in chromosomal_sequences:
            chromosomal_basecount += sequence.get_basecount()
        for sequence in extrachromosomal_sequences:
            extrachromosomal_basecount += sequence.get_basecount()

        return chromosomal_basecount / extrachromosomal_basecount

    
    def get_chromosome_sequences(self):
        return [x for x in self.sequences if 'plasmid' not in x.name and 'mitochond' not in x.name]
    

    def get_extrachromosomal_sequences(self):
        return [x for x in self.sequences if 'plasmid' in x.name or 'mitochond' in x.name]


    def get_hmq_reads(self):
        mapq_60_read_ids = []
        mapq_10_read_ids = []
        for seq in self.sequences:
            mapq_60_read_ids += seq.get_mapq_read_ids(60)
            mapq_10_read_ids += seq.get_mapq_read_ids(10)
        return mapq_60_read_ids, mapq_10_read_ids



class Sequence:
    def __init__(self, accession, name):
        self.accession = accession
        self.name = name
        self.alignments = [] # [read id, read length, base matches, block length, mapq, bases awarded]
    

    def get_basecount(self):
        if not hasattr(self, 'basecount'):
            alignments = self.alignments
            self.basecount = sum([x[5] for x in alignments])
        return self.basecount


    def get_read_ids(self):
        if not hasattr(self, 'read_ids'):
            alignments = self.alignments
            self.read_ids = [x[6] for x in alignments]
        return self.read_ids


    def get_hma_count(self):
        if not hasattr(self, 'mapq_60_count'):
            alignments = self.alignments
            self.mapq_60_count = sum([x[4] >= 60 for x in alignments])
        if not hasattr(self, 'mapq_10_count'):
            alignments = self.alignments
            self.mapq_10_count = sum([x[4] >= 10 for x in alignments])
        return self.mapq_60_count, self.mapq_10_count

    
    def get_mapq_read_ids(self, mapq_score):
        alignments = self.alignments
        return [x[9] for x in alignments if x[4] >= mapq_score]