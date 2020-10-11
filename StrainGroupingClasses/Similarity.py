
from collections import defaultdict



class SimilarityCalculator:
    def __init__(self, read_classifications):
        self.read_classifications = read_classifications


    def pairwise_occurance(self):
        self.get_raw_co_occurances()
        self.format_as_rate()
        return self.read_counts, self.co_occurance


    def get_raw_co_occurances(self):
        read_classifications = self.read_classifications
        co_occurance = {}
        read_counts = defaultdict(int)

        for read in read_classifications:
            identifiers = read.classifications 
            for ident in identifiers:
                read_counts[ident] += 1
                if ident not in co_occurance:
                    co_occurance[ident] = defaultdict(int)

            for i in range(len(identifiers)):
                for j in range(i + 1, len(identifiers)):
                    co_occurance[identifiers[i]][identifiers[j]] += 1
                    co_occurance[identifiers[j]][identifiers[i]] += 1

        self.read_counts = read_counts 
        self.co_occurance = co_occurance 
        

    def format_as_rate(self):
        read_counts = self.read_counts
        co_occurance = self.co_occurance

        for identifier in co_occurance:
            total_classifications = read_counts[identifier]
            for partner, count in co_occurance[identifier].items():
                co_occurance[identifier][partner] = count / total_classifications
        
        self.co_occurance = co_occurance








