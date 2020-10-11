
from modules.utils import load_banlist
import os


class ReferenceMaker:
    def __init__(self, strain_groups, context):
        self.strain_groups = strain_groups
        self.project_path = context.project_path
        self.database_path = context.database_path
        self.banlist = load_banlist(self.project_path + '/banlist.txt')
        self.min_strain_basecount = 1 * 100000


    def make(self):
        print('creating strain group reference metagenomes')
        for group in self.strain_groups:
            relevant_filenames = self.get_relevant_genome_files(group)
            self.write_database(relevant_filenames, group.id)


    def get_relevant_genome_files(self, group):
        relevant_filenames = []
        for filename, basecount in group.strain_basecounts.items():
            if basecount > self.min_strain_basecount:
                relevant_filenames.append(filename)

        relevant_filenames = [x for x in relevant_filenames if x not in self.banlist]
        relevant_filenames = [self.database_path + filename for filename in relevant_filenames]

        return relevant_filenames


    def write_database(self, relevant_filenames, group_id):
        folder = self.project_path + '/runtimefiles/group_databases/'
        database_filename = folder + 'strain_group_' + str(group_id) + '_initial.fa'
        
        with open(database_filename, 'w') as fp:
            for genome_path in relevant_filenames:
                if not os.path.exists(genome_path):
                    print(f'could not locate {genome_path}')
                    continue
                with open(genome_path, 'r') as fp2:
                    data = fp2.read()
                    fp.write(data)

