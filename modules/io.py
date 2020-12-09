
import json
from AbundanceClasses.Genome import Genome
import os
from collections import defaultdict

 

def load_json(filename):
    with open(filename, 'r') as fp:
        return json.load(fp)


def remove_reads_from_fastq(reads_filename, read_ids):
    if len(read_ids) < 10:
        print('no human reads')
        return

    read_ids = ['@' + r for r in read_ids]
    read_ids = set(read_ids)
    with open(reads_filename, 'r') as fq_pointer:
        newfile = reads_filename[:-3] + '_no_human.fq'
        with open(newfile, 'w') as new_fq_pointer:
            while 1:
                try:
                    line = fq_pointer.readline()
                    delim_loc = line.find(' ')
                    identifier = line[:delim_loc]

                    if identifier not in read_ids:
                        read = [line] + [next(fq_pointer) for i in range(3)]
                        new_fq_pointer.writelines(read)
                    
                    else:
                        next(fq_pointer)
                        next(fq_pointer)
                        next(fq_pointer)
                except StopIteration:
                    return


def load_genomes(folder):
    genomes = []
    fasta_files = get_fastas_in_folder(folder)
    fasta_files = [folder + f for f in fasta_files]
    for f in fasta_files:
        new_genome = load_genome(f)
        genomes.append(new_genome)
    return genomes


def load_genome(filename):
    genome = Genome(filename)
    with open(filename, 'r') as fp:
        active_header = ''
        active_sequence = ''
        line = fp.readline()
        while line:
            if line.startswith('>'):
                if active_header != '' and active_sequence != '':
                    genome.add_sequence(active_header, active_sequence)
                active_header = line.rstrip('\n')
                active_sequence = ''
            else:
                active_sequence += line.rstrip('\n')
            line = fp.readline()
        genome.add_sequence(active_header, active_sequence)
    return genome


def get_fastas_in_folder(folder):
    filenames = os.listdir(folder)
    fastas = []
    for f in filenames:
        extension = f.rsplit('.', 1)[1]
        if extension in ['fna', 'fasta', 'fa']:
            fastas.append(f)
    return fastas


def read_fasta(filename):
    with open(filename, 'r') as fp:
        sequence = ''
        header = fp.readline().rstrip('\n').lstrip('>')
        accession, genome_name = header.split(' ', 1)
        accession = accession.split('.', 1)[0]
        raw_sequence_list = fp.readlines()
        sequence_list = [line.rstrip('\n') for line in raw_sequence_list]
        sequence_list = [line for line in sequence_list if not line.startswith('>')]
        sequence = ''.join(sequence_list)
    return filename, accession, genome_name, sequence


def write_fastq(filenames_reads, filename):
    with open(filename, 'w') as fp:
        read_count = 0
        for filename, reads in filenames_reads.items():
            for read in reads:
                read_count += 1
                fp.write('@' + filename + '_' + str(read_count) + '\n')
                fp.write(read + '\n')
                fp.write('+\n')
                fp.write('-' * len(read) + '\n')


def write_json(the_dict, filename):
    with open(filename, 'w') as fp:
        json.dump(the_dict, fp)



def write_reads_from_dict(the_dict, fastq_file, output_filename):
    genomes_reads = []

    read_ids_genomes = {}
    for genome, read_ids in the_dict.items():
        for read_id in read_ids:
            read_ids_genomes[read_id] = genome

    with open(fastq_file, 'r') as read_fp:
        with open(output_filename, 'w') as out_fp:
            line = read_fp.readline()
            while line:
                if line.startswith('@'):
                    read_id = line.split(' ', 1)[0].lstrip('@').rstrip('\n')
                    if read_id in read_ids_genomes:
                        genome = read_ids_genomes[read_id]
                        sequence = read_fp.readline()
                        out_fp.write('\t'.join([genome, read_id, sequence]) + '\n')

                        line = read_fp.readline()
                        if line != '+':
                            continue
                        else:
                            line = read_fp.readline()

                line = read_fp.readline()



def pull_reads_by_id(read_ids, original_fastq, output_fastq):
    required_reads = set(read_ids)

    with open(original_fastq, 'r') as read_fp:
        with open(output_fastq, 'w') as out_fp:
            line = read_fp.readline()
            while line:
                if line.startswith('@'):
                    read_id = line.split(' ', 1)[0].lstrip('@').rstrip('\n')
                    if read_id in required_reads:
                        out_fp.write(line)
                        sequence = read_fp.readline()
                        plus = read_fp.readline()
                        qstring = read_fp.readline()
                        out_fp.write(sequence)
                        out_fp.write(plus)
                        out_fp.write(qstring)

                line = read_fp.readline()



def write_group_genomes(genome_list, filename):
    with open(filename, 'w') as fp:
        for genome in genome_list:
            fp.write(genome + '\n')

    
