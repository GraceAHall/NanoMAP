

from GeneralClasses.Clades import Strain, Sequence
from collections import defaultdict


class SampleCharacterisation:
    def __init__(self, af_dict, an_dict):
        self.name = 'sample'
        self.af_dict = af_dict
        self.an_dict = an_dict
        self.strains = {}
        self.multimapping_dict = {}
        self.sampled_base_count = 0
        self.min_cec_ratio = 2
        self.resolved = False


    def update_characterisation(self, read_id, alignments):
        af_dict = self.af_dict
        an_dict = self.an_dict
        read_length = alignments[0][1]
        equal_bases_share = read_length / len(alignments)
        strains = self.strains

        filename_classifications = set()
        for al in alignments:
            al.append(equal_bases_share)
            al.append(read_id)
            accession = al[0].upper()
            filename = af_dict[accession]
            sequencename = an_dict[accession].lower()
            filename_classifications.add(filename)
            
            if filename not in strains:
                strains[filename] = Strain(filename)
            
            if accession not in strains[filename].sequences:
                strains[filename].sequences[accession] = Sequence(accession, sequencename)
            
            strains[filename].sequences[accession].alignments.append(al)

        filename_classifications = list(filename_classifications)
        num_classifications = len(filename_classifications)
        if filename_classifications[0] not in self.multimapping_dict:
            self.multimapping_dict[filename_classifications[0]] = defaultdict(int) 

        for i in range(num_classifications - 1):
            f1 = filename_classifications[i]
            if f1 not in self.multimapping_dict:
                self.multimapping_dict[f1] = defaultdict(int)

            for j in range(i + 1, num_classifications):
                f2 = filename_classifications[j]
                if f2 not in self.multimapping_dict:
                    self.multimapping_dict[f2] = defaultdict(int)
                self.multimapping_dict[f1][f2] += 1
                self.multimapping_dict[f2][f1] += 1


    def reformat_strains_to_lists(self):
        self.strains = list(self.strains.values())
        for strain in self.strains:
            strain.sequences = list(strain.sequences.values())
            

    def remove_extrachromosomal_dominant_strains(self):
        strains_to_delete = []
        for i, strain in enumerate(self.strains):
            cec_ratio = strain.get_chromosomal_extrachromosomal_ratio()
            if cec_ratio < self.min_cec_ratio:
                strains_to_delete.append(i)
        
        print(f'{len(strains_to_delete)} strains are extrachromosomal dominant')
        for index in reversed(strains_to_delete):
            del self.strains[index]


    def set_sampled_base_count(self):
        base_count = 0
        assert(type(self.strains) == list)
        for strain in self.strains:
            base_count += strain.get_total_basecount()
        self.sampled_base_count = base_count
        

    def set_sample_abundances(self):
        for strain in self.strains:
            strain.group_abundance = strain.get_total_basecount() / self.sampled_base_count * 100


    def set_mapq_info(self):
        for strain in self.strains:
            strain.get_hma_counts() 

    
    def get_HMQ_read_ids(self):
        hmq60_reads = defaultdict(dict)
        hmq10_reads = defaultdict(dict)
        for strain in self.strains:
            name = strain.filename
            mapq_60_read_ids, mapq_10_read_ids = strain.get_hmq_reads()
            if len(mapq_10_read_ids) > 0:
                hmq60_reads[name] = mapq_60_read_ids
                hmq10_reads[name] = mapq_10_read_ids
        return hmq60_reads, hmq10_reads


    def print_strains(self):
        print(self.name)
        print('\n\n{:>100}{:>50}{:>10}{:>10}{:>10}{:>10}{:>10}{:>12}'.format('name', 'filename', 'naive', 'MAPQ_60', 'MAPQ_10', 'block ratio', 'presence', 'sample abund.'))
        self.strains.sort(key=lambda x: x.mapq_10_count, reverse=True)
        for s in self.strains:
            name = s.get_primary_name()
            name = ' '.join(name.split(' ')[-5:])
            filename = s.filename
            print('{:>100}{:>50}{:>10.2f}{:>10.1f}{:>10.1f}{:>10.1f}{:>10}{:>12.1f}'.format(name[:96], filename[:46], s.group_abundance, s.mapq_60_count, s.mapq_10_count, s.block_ratio, s.mapq_confirmed, s.sample_abundance))

    
