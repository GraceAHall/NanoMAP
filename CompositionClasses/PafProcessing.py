

from collections import defaultdict
from modules.io import load_json



class PafProcessor:
    def __init__(self, paf_file, database_path):
        self.paf_file = paf_file
        self.ag_dict = load_json(database_path + 'taxonomy/accessions_genomes.json')
        self.af_dict = load_json(database_path + 'taxonomy/accessions_filenames.json')


    def process(self):
        # [read_id, target, collinearity, pid, block, read_length, base_matches, MAPQ, start_offset, end_offset]
        alignments = self.load_alignments()
        reads_alignments = self.group_by_read_id(alignments)
        reads_alignments = self.pick_best_alignments(reads_alignments)
        return list(reads_alignments.values())


    def load_alignments(self):
        with open(self.paf_file, 'r') as fp:
            lines = fp.readlines()

        alignments = []
        for line in lines:
            line = line.split('\t')
            if len(line) > 1:
                read_id = line[0]
                target = line[5]
                read_length = int(line[1])
                query_len = int(line[3]) - int(line[2])
                target_len = int(line[8]) - int(line[7])
                collinearity = 1 - abs((query_len - target_len) / query_len)
                pid = int(line[9]) / int(line[10])
                block = int(line[9])
                base_matches = int(line[9])
                mapq = int(line[11])
                start_offset = int(line[7])
                end_offset = int(line[6]) - int(line[8])
                alignments.append([read_id, target, collinearity, pid, block, read_length, base_matches, mapq, start_offset, end_offset])

        del lines
        return alignments


    def get_collinearities(self, alignments):
        collinearities = defaultdict(list)
        af_dict = self.af_dict

        for al in alignments:
            filename = af_dict[al[1].upper()]
            collinearities[filename].append(al[2])

        for filename, collinearity_list in collinearities.items():
            collinearities[filename] = sum(collinearity_list) / len(collinearity_list) * 100

        return collinearities


    def group_by_read_id(self, alignments):
        reads_alignments = defaultdict(list)
        for al in alignments:
            reads_alignments[al[0]].append(al)
        return reads_alignments

    
    def pick_best_alignments(self, reads_alignments):
        for read_id, alignments in reads_alignments.items():
            if len(alignments) != 0:
                read_length = alignments[0][5]
                highest_mapq = max([x[7] for x in alignments])
                alignments = [x for x in alignments if x[7] == highest_mapq]

                alignments.sort(key=lambda x: x[6], reverse=True)
                reads_alignments[read_id] = alignments[0]
        
        return reads_alignments