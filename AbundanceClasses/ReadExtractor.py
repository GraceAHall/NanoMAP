

import math
import numpy as np
import os
from subprocess import Popen


class ReadExtractor:
    def __init__(self, genomes):
        self.genomes = genomes
        self.read_length = 3000


    def extract(self, method='stride', accessory_dna=False):
        reads = {}
        if method == 'stride':
            for genome in self.genomes:
                reads[genome.accession] = self.get_reads_stride(genome)
        elif method == 'randunif':
            for genome in self.genomes:
                reads[genome.accession] = self.get_reads_randunif(genome)
        elif method == 'badread':
            self.get_reads_badread(self.genomes)
        return reads


    def get_reads_stride(self, genome):
        reads = []
        num_reads = math.floor(genome.basecount / self.read_length)
        stride = math.ceil(genome.length / num_reads)
        sequence = genome.sequence + genome.sequence[:100000]
        for i in range(num_reads):
            extract_location = i * stride
            if extract_location > len(sequence):
                break
            reads.append(sequence[extract_location: extract_location + self.read_length])
        return reads


    def get_reads_randunif(self, genome):
        reads = []
        num_reads = math.floor(genome.basecount / self.read_length)
        sequence = genome.sequence + genome.sequence[:3000]
        extract_locations = np.random.uniform(0, len(sequence) - 3000, num_reads)
        for loc in extract_locations:
            loc = int(loc)
            reads.append(sequence[loc: loc + self.read_length])
        return reads
    

    def get_reads_badread(self, genomes):
        print('simulating reads with badread')
        if os.path.exists('programfiles/reads.fq'):
            os.remove('programfiles/reads.fq')

        processes = []
        for i in range(len(genomes)):
            million_bases = genomes[i].basecount // 1000000
            print(million_bases)
            with open(f'programfiles/reads_{i}.fq', 'w') as outfile:
                cmd = f'../Badread/badread-runner.py simulate --reference {genomes[i].filename} --quantity {million_bases}M \
                    --length 3000,0 --identity 88,95,5 --start_adapter 0,0 --end_adapter 0,0 --junk_reads \
                    0.5 --random_reads 0.5 --chimeras 0.5'
                p = Popen(cmd, shell=True, stdout=outfile) 
                processes.append(p)
        
        for p in processes:
            p.wait()