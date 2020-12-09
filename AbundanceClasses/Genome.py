
 
class Genome:
    def __init__(self, filename):
        self.kmersize = 15
        self.filename = filename.rsplit('/', 1)[1]
        self.sequences = []


    def add_sequence(self, header, sequence):
        new_sequence = Sequence(header, sequence)
        self.sequences.append(new_sequence)


    def get_accessions(self):
        accessions = []
        for sequence in self.sequences:
            accessions.append(sequence.accession)
        return accessions


    def get_chromsome_sequences(self):
        chromosome_sequences = []
        for sequence in self.sequences:
            if 'plasmid' not in sequence.header.lower() and 'mitochon' not in sequence.header.lower():
                chromosome_sequences.append(sequence)
        return chromosome_sequences


    def get_all_dna(self):
        dna = ''
        for seq in self.sequences:
            dna += seq.sequence
        return dna

    
    def get_longest_sequence(self):
        self.sequences.sort(key=lambda x: x.length, reverse=True)
        return self.sequences[0]
            

    @staticmethod
    def generate_reverse_compliment(sequence):    
        reverse_sequence = sequence[::-1]
        origin = 'ACGT'
        desination = 'TGCA'
        table = reverse_sequence.maketrans(origin, desination)
        revseq = reverse_sequence.translate(table) 
        return revseq


class Sequence:
    def __init__(self, header, sequence):
        self.header = header
        self.accession = header.split(' ', 1)[0].lstrip('>')
        self.sequence = sequence
        self.length = len(sequence)
