



from collections import defaultdict
from modules.io import load_json



class PafProcessor:
    def __init__(self, paf_file, database_path, include_plasmids_mitochondria):
        self.paf_file = paf_file
        self.include_plasmids_mitochondria = include_plasmids_mitochondria
        self.ag_dict = load_json(database_path + 'taxonomy/accessions_genomes.json')
        self.af_dict = load_json(database_path + 'taxonomy/accessions_filenames.json')


    def process(self):
        # [read_id, target, collinearity, pid, block, read_length, base_matches, MAPQ, start_offset, end_offset]
        alignments = self.load_alignments()
        collinearities = self.get_collinearities(alignments)
        reads_alignments = self.group_by_read_id(alignments)
        reads_alignments = self.remove_origin_mapping_reads(reads_alignments)
        reads_alignments = self.filter_length(reads_alignments)
        reads_alignments = self.filter_short_alignments(reads_alignments)
        reads_alignments = self.filter_pid(reads_alignments)
        reads_alignments = self.filter_non_collinear(reads_alignments)
        if not self.include_plasmids_mitochondria:
            reads_alignments = self.remove_extrachromosomal_alignments(reads_alignments)
        reads_alignments = self.pick_best_alignments(reads_alignments)
        return reads_alignments, collinearities


    def load_alignments(self):
        with open(self.paf_file, 'r') as fp:
            lines = fp.readlines()
            lines = [ln.split('\t') for ln in lines]
            lines = [ln for ln in lines if len(ln) > 1]

        alignments = []
        for line in lines:
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

        avg_collinearity = sum([al[2] for al in alignments]) / len(alignments)
        avg_pid = sum([al[3] for al in alignments]) / len(alignments)

        std_dev_collinearity = sum([(al[2] - avg_collinearity) ** 2 for al in alignments]) / (len(alignments) - 1)
        std_dev_pid = sum([(al[3] - avg_pid) ** 2 for al in alignments]) / (len(alignments) - 1)

        self.min_collinearity = avg_collinearity - 1 * std_dev_collinearity
        self.min_pid = avg_pid - 1 * std_dev_pid

        #print('min collinearity', self.min_collinearity)
        #print('min pid', self.min_pid)
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


    def remove_origin_mapping_reads(self, reads_alignments):
        reads_to_delete = set()
        for read_id, alignments in reads_alignments.items():
            for al in alignments:
                if al[8] < 100 or al[9] < 100:
                    reads_to_delete.add(read_id)

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]
        
        return reads_alignments


    def filter_length(self, reads_alignments):
        min_length = 1000
        reads_to_delete = set()

        for read_id, alignments in reads_alignments.items():
            if alignments[0][5] < min_length:
                reads_to_delete.add(read_id)

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]

        return reads_alignments


    def filter_short_alignments(self, reads_alignments):
        min_length = 1000
        reads_to_delete = set()

        for read_id, alignments in reads_alignments.items():
            filtered_alignments = [al for al in alignments if al[4] > min_length]
            if len(filtered_alignments) == 0:
                reads_to_delete.add(read_id)
            else:
                reads_alignments[read_id] = filtered_alignments

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]

        return reads_alignments


    def filter_pid(self, reads_alignments):
        min_pid = self.min_pid
        reads_to_delete = set()

        for read_id, alignments in reads_alignments.items():
            filtered_alignments = [al for al in alignments if al[3] > min_pid]
            if len(filtered_alignments) == 0:
                reads_to_delete.add(read_id)
            else:
                reads_alignments[read_id] = filtered_alignments

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]

        return reads_alignments
            

    def filter_non_collinear(self, reads_alignments):
        min_collinearity = self.min_collinearity
        reads_to_delete = set()

        for read_id, alignments in reads_alignments.items():
            filtered_alignments = [al for al in alignments if al[2] > min_collinearity]
            if len(filtered_alignments) == 0:
                reads_to_delete.add(read_id)
            else:
                reads_alignments[read_id] = filtered_alignments

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]

        return reads_alignments

    
    def remove_extrachromosomal_alignments(self, reads_alignments):
        #print('\n\nremoving extrachromosomal\n\n')
        reads_to_delete = set()
        taxonomy = self.ag_dict

        for read_id, alignments in reads_alignments.items():
            filtered_alignments = []
            for al in alignments:
                target = taxonomy[al[1].upper()]
                if 'plasmid' not in target and 'mitochond' not in target:
                    filtered_alignments.append(al)

            if len(filtered_alignments) == 0:
                reads_to_delete.add(read_id)
            else:
                reads_alignments[read_id] = filtered_alignments

        for read_id in list(reads_to_delete):
            del reads_alignments[read_id]

        return reads_alignments


    def pick_best_alignments(self, reads_alignments):
        for read_id, alignments in reads_alignments.items():
            if len(alignments) != 0:
                read_length = alignments[0][5]
                highest_mapq = max([x[7] for x in alignments])
                best_mapqs = [x for x in alignments if x[7] == highest_mapq]
                highest_base_matches = max([x[6] for x in best_mapqs])
                best_mapqs_good_matches = [x for x in best_mapqs if x[6] > 0.999 * highest_base_matches]
                reads_alignments[read_id] = best_mapqs_good_matches
        
        return reads_alignments





