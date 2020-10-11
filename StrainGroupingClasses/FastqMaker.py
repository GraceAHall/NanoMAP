

from collections import defaultdict


class FastqMaker:
    def __init__(self, strain_groups, context):
        self.strain_groups = strain_groups
        self.project_path = context.project_path
        self.fastq_path = context.fastq_path


    def make(self):
        self.set_file_mode()
        self.make_read_id_dict()

        if self.mode == 'fasta':
            print('binning reads from input fasta')
            self.set_reads_by_group_fasta()
        else:
            print('binning reads from input fastq')
            self.set_reads_by_group_fastq()

        self.write_new_fastqs()
 
    
    def set_file_mode(self):
        with open(self.fastq_path, 'r') as fp:
            line = fp.readline()
            header_mark = line[0]
            mode = 'fasta' if header_mark == '>' else 'fastq'
        self.mode = mode


    def make_read_id_dict(self):
        banlist = set()
        read_id_dict = {}
        for group in self.strain_groups:
            for read_id in list(group.read_ids):
                if read_id in read_id_dict:
                    del read_id_dict[read_id]
                    banlist.add(read_id)
                if read_id not in banlist:
                    read_id_dict[read_id] = group.id
        print(f'{len(banlist)} multigroup reads removed')
        self.read_id_dict = read_id_dict


    def set_reads_by_group_fastq(self):
        read_id_dict = self.read_id_dict
        group_reads = defaultdict(list)

        with open(self.fastq_path, 'r') as fp:
            line = fp.readline()
            while line:
                if line.startswith('@'):
                    read_id = line.split(' ')[0].lstrip('@').rstrip('\n')
                    if read_id in read_id_dict:
                        group = read_id_dict[read_id]
                        group_reads[group].append(line)
                        group_reads[group].append(fp.readline())
                        group_reads[group].append(fp.readline())
                        group_reads[group].append(fp.readline())
                    else:
                        next(fp)
                        next(fp)
                        next(fp)
                line = fp.readline()
        self.group_reads = group_reads


    def set_reads_by_group_fasta(self):
        read_id_dict = self.read_id_dict
        group_reads = defaultdict(list)

        header = ''
        sequence = ''

        with open(self.fastq_path, 'r') as fp:
            line = fp.readline()
            while line:
                if line.startswith('>'):
                    if header != '':
                        read_id = header.split(' ')[0].lstrip('>').rstrip('\n')
                        if read_id in read_id_dict:
                            group = read_id_dict[read_id]
                            group_reads[group].append('@' + header[1:])
                            group_reads[group].append(sequence + '\n')
                            group_reads[group].append('+\n')
                            group_reads[group].append('-' * len(sequence) + '\n')

                    header = line
                    sequence = ''
                else:
                    sequence += line.rstrip('\n')

                line = fp.readline()

        read_id = header.split(' ')[0].lstrip('>').rstrip('\n')
        if read_id in read_id_dict:
            group = read_id_dict[read_id]
            group_reads[group].append('@' + header[1:])
            group_reads[group].append(sequence + '\n')
            group_reads[group].append('+\n')
            group_reads[group].append('-' * len(sequence) + '\n')

        self.group_reads = group_reads


    def write_new_fastqs(self):
        group_reads = self.group_reads
        folder = self.project_path + '/runtimefiles/group_fastqs/'

        for group_id, lines in group_reads.items():
            group_filename = folder + 'strain_group_' + str(group_id) + '.fq'
            with open(group_filename, 'w') as fp:
                fp.writelines(lines) 


