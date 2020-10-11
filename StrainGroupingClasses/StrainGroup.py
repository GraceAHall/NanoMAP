


from collections import defaultdict



class StrainGroup:
    def __init__(self, id):
        self.id = id
        self.basecount = 0
        self.abundance = 0
        self.read_ids = set()
        self.identifiers = set()
        self.strain_basecounts = defaultdict(int)