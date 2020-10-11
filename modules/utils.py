
#import matplotlib.pyplot as plt
#import numpy as np



def load_banlist(filename):
    with open(filename, 'r') as fp:
        lines = fp.readlines()
        lines = [x.rstrip('\n') for x in lines]
        banlist = set(lines)
        return banlist



def plot_fastq_read_lengths(fastq_filename):
    readlengths = []
    with open(fastq_filename, 'r') as fp:
        line = fp.readline()
        while line:
            if line.startswith('@'):
                sequence = fp.readline()
                readlengths.append(len(sequence))
                next(fp)
                next(fp)
            line = fp.readline()

    print(f'avg read length: {sum(readlengths) / len(readlengths)}')
    bins = np.arange(0, 20000, 100)
    plt.hist(readlengths, bins=bins)
    plt.title(f'read length distribution in {fastq_filename}\ntotal reads: {len(readlengths)}')
    plt.xlabel('read length')
    plt.ylabel('count')
    plt.show()




def get_best_alignment(lines):
    best_alignments = []
    best = lines[0].split('\t')
    best[1] = int(best[1])
    best[9] = int(best[9])
    best[10] = int(best[10])

    for line in lines:
        details = line.split('\t')
        details[1] = int(details[1])
        details[9] = int(details[9])
        details[10] = int(details[10])

        if details[0] == best[0]:
            if details[9] > best[9]:
                #print('\nnew best:')
                #print(f'old: {best[:12]}')
                #print(f'new: {details[:12]}')
                best = details
        else:
            best_alignments.append([best[5], int(best[1]), int(best[9]), int(best[10]), int(best[11])])
            best = details
            best[1] = int(best[1])
            best[9] = int(best[9])
            best[10] = int(best[10])
            
    best_alignments.append([best[5], int(best[1]), int(best[9]), int(best[10]), int(best[11])])
    return best_alignments








def filter_low_mapq_reads(alignments, min_mapq):
    # accession, read len, num residue matches, block length, mapq
    output = []
    for al in alignments:
        if al[4] >= min_mapq:
            output.append(al)

    return output


def filter_short_reads(alignments, min_read_length):
    # accession, read len, num residue matches, block length, mapq
    output = []
    for al in alignments:
        if al[1] >= min_read_length:
            output.append(al)

    return output






def update_count_dict(alignments, the_dict):
    for al in alignments:
        the_dict[al[0]] += 1
    return the_dict


def map_accession_to_genome(accessions_counts, taxonomy):
    accessions_counts = list(accessions_counts.items())
    
    translated = {}
    for item, count in accessions_counts:
        try:
            translated[taxonomy[item.upper()]] = count
        except KeyError:
            translated[item.upper()] = count

    return translated


def group_by_species(abundances):
    output = {}

    for genome, count in abundances.items():
        species = genome.split(' ', 2)[:2]
        species = ' '.join(species).lower()
        species = species.rstrip(',')
        if species not in output:
            output[species] = []
        
        output[species].append([genome, count])

    return output


def filter_short_alignments(alignments, min_alignment_percent):
    # accession, read len, num residue matches, block length, mapq
    passing_alignments = []
    for al in alignments:
        if al[4] > al[1] * min_alignment_percent:
            passing_alignments.append(al)

    #print(f'start len: {len(alignments)} reads')
    #print(f'end len: {len(passing_alignments)} reads')
    return passing_alignments


def filter_low_identity_alignments(alignments, min_pid):
    # accession, read len, num residue matches, block length, mapq, primary/secondary
    passing_alignments = []
    for al in alignments:
        if al[3] > al[4] * min_pid:
            passing_alignments.append(al)
            
    return passing_alignments

