

from collections import defaultdict
from modules.io import load_json
from StrainGroupingClasses.ReadClassification import ReadClassification



class PafProcessor:
    def __init__(self, paf_file, database_path):
        self.paf_file = paf_file
        self.min_pid = 0.1
        self.min_collinearity = 0.8
        self.min_read_length = 700
        self.ag_dict = load_json(database_path + 'taxonomy/accessions_genomes.json')
        self.af_dict = load_json(database_path + 'taxonomy/accessions_filenames.json')


    def process(self):
        # [read_id, target, collinearity, pid, block, read_length, base_matches]
        alignments = self.load_alignments()
        alignments = self.filter_read_length(alignments)
        alignments = self.filter_pid(alignments)
        alignments = self.filter_non_collinear(alignments)
        reads_alignments = self.group_by_read_id(alignments)
        reads_alignments = self.translate_to_filenames(reads_alignments)
        strains_basecounts = self.get_strain_basecounts(reads_alignments)
        read_classifications = self.classify_reads(reads_alignments)
        return read_classifications, strains_basecounts


    def load_alignments(self):
        alignments = []
        with open(self.paf_file, 'r') as fp:
            line = fp.readline().split('\t')
            while line and len(line) > 1:
                read_id = line[0]
                target = line[5]
                read_length = int(line[1])
                query_len = int(line[3]) - int(line[2])
                target_len = int(line[8]) - int(line[7])
                collinearity = 1 - abs((query_len - target_len) / query_len)
                pid = int(line[9]) / int(line[10])
                block = int(line[9])
                base_matches = int(line[9])
                alignments.append([read_id, target, collinearity, pid, block, read_length, base_matches])
                line = fp.readline().split('\t')

        #print('collinearity', sum([al[2] for al in alignments]) / len(alignments))
        #print('pid', sum([al[3] for al in alignments]) / len(alignments))
        return alignments        



    def filter_read_length(self, alignments):
        min_read_length = self.min_read_length
        start_len = len(alignments)
        alignments = [al for al in alignments if al[5] > min_read_length]
        end_len = len(alignments)
        #print(f'removed {start_len - end_len} alignments due to read length')
        return alignments


    def filter_pid(self, alignments):
        min_pid = self.min_pid
        start_len = len(alignments)
        alignments = [al for al in alignments if al[3] > min_pid]
        end_len = len(alignments)
        #print(f'removed {start_len - end_len} alignments due to low pid')
        return alignments
            

    def filter_non_collinear(self, alignments):
        min_collinearity = self.min_collinearity
        start_len = len(alignments)
        alignments = [al for al in alignments if al[2] > min_collinearity]
        end_len = len(alignments)
        #print(f'removed {start_len - end_len} alignments due to low collinearity')
        return alignments 

    
    def group_by_read_id(self, alignments):
        reads_alignments = defaultdict(list)
        for al in alignments:
            reads_alignments[al[0]].append(al)
        return reads_alignments


    def translate_to_filenames(self, reads_alignments):
        for read_id, alignments in reads_alignments.items():
            if len(alignments) != 0:
                for i in range(len(alignments)):
                    try:
                        alignments[i][1] = self.af_dict[alignments[i][1].upper()]
                    except KeyError:
                        print(f'taxonomy error: {alignments[i][1]}')
        return reads_alignments
        

    def get_strain_basecounts(self, reads_alignments):
        strains_basecounts = defaultdict(int)
        for read_id, alignments in reads_alignments.items():
            if len(alignments) != 0:
                max_base_matches = max([x[6] for x in alignments])
                classifications = [x[1] for x in alignments if x[6] == max_base_matches]
                classifications = list(set(classifications))

                read_length = alignments[0][5]
                bases_share = read_length / len(classifications)

                for identifier in classifications:
                    strains_basecounts[identifier] += bases_share
                
        return strains_basecounts


    def classify_reads(self, reads_alignments):
        read_classifications = []
        for read_id, alignments in reads_alignments.items():
            if len(alignments) != 0:
                read_length = alignments[0][5]
                max_block = max([x[4] for x in alignments])
                good_block = [x for x in alignments if x[4] > 0.5 * max_block]
                max_pid = max([x[3] for x in good_block])
                good_block_good_pid = [x for x in good_block if x[3] > 0.9 * max_pid]
                identifiers = list(set([x[1] for x in good_block_good_pid]))
                new_classification = ReadClassification(read_id, read_length, identifiers)
                read_classifications.append(new_classification)
        return read_classifications



