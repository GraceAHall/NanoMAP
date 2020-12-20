

from collections import defaultdict
from modules.io import load_json, write_json
import itertools


class PafClassifier:
    def __init__(self, paf_file, accessions_filenames, symbol_mappings, context):
        self.paf_file = paf_file
        self.accessions_filenames = accessions_filenames
        self.symbol_mappings = symbol_mappings
        self.context = context


    def classify(self, empirical=False):
        self.group_alignments_by_read()
        self.classify_reads()
        self.remove_plasmids_mitochondria()
        self.translate_accessions_to_symbols()
        self.set_classification_combinations()
        if empirical:
            self.translate_read_ids_to_symbols()
            self.summarise_empirical()
            self.jsonify_empirical_counts()
        else:
            self.summarise_observed()
            self.jsonify_observed_counts()

        #self.print_classification_stats()


    def jsonify_empirical_counts(self):
        path = self.context.project_path + '/runtimefiles/abundance_estimation/'
        write_json(self.organisms_summary_dict, path + 'empirical_counts.json')
        write_json(self.pooled_summary_dict, path + 'empirical_pooled.json')


    def jsonify_observed_counts(self):
        path = self.context.project_path + '/runtimefiles/abundance_estimation/'
        write_json(self.pooled_summary_dict, path + 'observed_pooled.json')


    def group_alignments_by_read(self):
        filename = self.paf_file
        reads_alignments = defaultdict(list)
        with open(filename, 'r') as fp:
            line = fp.readline()
            while line:
                al = line.split('\t')
                reads_alignments[al[0]].append([al[5], int(al[1]), int(al[9]), int(al[10])])
                line = fp.readline()
        self.reads_alignments = reads_alignments


    def classify_reads(self):
        classifications = []
        for read_id, alignments in self.reads_alignments.items():
            alignments = self.get_good_blocks(alignments)
            alignments = self.filter_low_pid(alignments)
            alignments = self.get_best_base_matches(alignments)
            read_classification = self.format_to_final_alignments(alignments, read_id)
            if len(alignments) > 0:
                classifications.append(read_classification)
        self.classifications = classifications


    @staticmethod
    def get_good_blocks(alignments):
        if len(alignments) == 0:
            return alignments

        read_length = alignments[0][1]
        max_block = 1.1 * read_length + 50   # - constant for adaptor size?
        min_block = 0.9 * read_length - 50   # + constant for adaptor size?
        return [x for x in alignments if x[3] > min_block and x[3] < max_block]


    @staticmethod
    def filter_low_pid(alignments):
        if len(alignments) == 0:
            return alignments

        min_identity = 0.85
        return [al for al in alignments if al[2] / al[3] > min_identity]


    @staticmethod
    def get_best_base_matches(alignments):
        if len(alignments) == 0:
            return alignments

        alignments.sort(key=lambda x: x[2], reverse=True)  # sort by base matches
        highest_num_matches = alignments[0][2]
        return [x for x in alignments if x[2] == highest_num_matches]


    @staticmethod
    def format_to_final_alignments(alignments, read_id):
        if len(alignments) == 0:
            return alignments
        
        classifications = set()
        for al in alignments:
            classifications.add(al[0])
        classifications = [read_id, list(classifications)]
        return classifications


    def remove_plasmids_mitochondria(self):
        ct = self.context
        taxonomy_file = ct.database_path + 'taxonomy/accessions_genomes.json' 
        taxonomy = load_json(taxonomy_file)

        clean_classifications = []
        while len(self.classifications) > 0:
            chromosome_accessions = []
            read_id, accessions = self.classifications.pop()

            for accession in accessions:
                sequence_name = taxonomy[accession.upper()]
                if 'plasmid' not in sequence_name.lower() and 'mitochond' not in sequence_name.lower():
                    chromosome_accessions.append(accession)

            if len(chromosome_accessions) > 0:
                clean_classifications.append([read_id, chromosome_accessions])
        self.classifications = clean_classifications


    def translate_read_ids_to_symbols(self):
        output = []
        while len(self.classifications) > 0:
            read_id, classifications = self.classifications.pop()
            filename = read_id.rsplit('_', 1)[0]
            symbol = self.symbol_mappings[filename]
            output.append([symbol, classifications])
        self.classifications = output

    
    def translate_accessions_to_symbols(self):
        output = []
        while len(self.classifications) > 0:
            identifier, classifications = self.classifications.pop()
            symbol_classifications = []
            for accession in classifications:
                filename = self.accessions_filenames[accession]
                symbol = self.symbol_mappings[filename]
                symbol_classifications.append(symbol)
            output.append([identifier, symbol_classifications])
        self.classifications = output


    def summarise_empirical(self):
        organisms_summary_dict = self.initialise_organism_summary_dict()
        self.create_organism_summaries(organisms_summary_dict)
        self.create_pooled_summary()


    def initialise_organism_summary_dict(self):
        symbols = list(self.symbol_mappings.values())

        organisms_summary_dict = {}
        for s in symbols:
            organisms_summary_dict[s] = {}
            for item in self.combinations:
                organisms_summary_dict[s][item] = 0

        return organisms_summary_dict


    def set_classification_combinations(self):
        symbols = list(self.symbol_mappings.values())

        # check non-empty list
        if len(symbols) == 0:
            self.combinations = []
            return

        classification_combinations = []
        for i in range(1, len(symbols) + 1):
            combinations = list(itertools.combinations(symbols, i))
            for comb in combinations:
                label = ''.join(comb)
                classification_combinations.append(label)
        self.combinations = classification_combinations


    def create_organism_summaries(self, organisms_summary_dict): 
        for organism, classifications in self.classifications: 
            classifications = list(set(classifications))
            code = ''.join(sorted(classifications))
            organisms_summary_dict[organism][code] += 1
        self.organisms_summary_dict = organisms_summary_dict
        

    def create_pooled_summary(self):
        pooled_summary = {}

        for item in self.combinations:
            pooled_summary[item] = 0

        for organism, summary in self.organisms_summary_dict.items():
            for classification, count in summary.items():
                pooled_summary[classification] += count

        self.pooled_summary_dict = pooled_summary   

    
    def summarise_observed(self):
        symbols_counts = self.initialise_symbols_counts()
        self.create_summary(symbols_counts)


    def initialise_symbols_counts(self):
        symbols_counts = {}
        for c in self.combinations:
            symbols_counts[c] = 0
 
        return symbols_counts


    def create_summary(self, symbols_counts):
        for read_id, classifications in self.classifications: 
            classifications = list(set(classifications))
            code = ''.join(sorted(classifications))
            symbols_counts[code] += 1
        self.pooled_summary_dict = symbols_counts


    def print_classification_stats(self):
        print('\n\n########## Pooled summary ##########')
        for symbol, count in self.pooled_summary_dict.items():
            print(f'classification: {symbol} count: {count}')
        
        if hasattr(self, 'organisms_summary_dict'):
            print('\n\n########## Organism summary ##########')
            for organism, symbol_counts in self.organisms_summary_dict.items():
                for symbol, count in symbol_counts.items():
                    print(f'origin: {organism} classification: {symbol} count: {count}')


