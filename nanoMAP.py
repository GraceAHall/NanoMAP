

from GeneralClasses.Context import ProgramContext, StrainGroupingContext, GroupContext
from GeneralClasses.FileSetup import FileSetup
from GeneralClasses.Minimap2Runner import Minimap2Runner
from StrainGroupingClasses.StrainGrouping import StrainGrouper
from GroupClasses.GroupAnalysis import GroupAnalyser
from CompositionClasses.CompositionAnalysis import CompositionAnalyser
#from GeneralClasses.Reporter import Reporter

import sys
import getopt
import subprocess
import os



def main(argv):
    opts, args = getopt.getopt(argv, "ht:r:d:p:m:", ["map-pb", "map-ont", "no-initial-alignment", "allow-extrachromosomal"])
    context = ProgramContext(opts) 
    check_environment()
    check_database_build(context)

    params = [opt[0] for opt in opts]
    if "--no-initial-alignment" in params:
        perform_noalign_setup(context)
    else:
        perform_full_setup(context)
        initial_alignment(context)

    strain_groups = group_strains(context)
    identified_strains = identify_sample_strains(strain_groups, context)
    report_results(identified_strains, context)
                           

def check_environment():
    if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 6):
        raise Exception("Must be using Python 3.6 or greater")


def check_database_build(context):
    try:
        path = os.path.exists(context.database_path + 'taxonomy/accessions_filenames.json')
    except FileNotFoundError:
        print('Could not load database. Has it been built?')
        sys.exit()


def perform_full_setup(context):
    fs = FileSetup(context)
    fs.full_setup()


def perform_noalign_setup(context):
    fs = FileSetup(context)
    fs.noalign_setup()


def initial_alignment(context):
    print('\n\n=========== Full Database Alignment ===========')
    mm2 = Minimap2Runner(context)
    mm2.initial_alignment()
    

def group_strains(context):
    print('\n\n=========== Strain Grouping ===========')
    grouping_context = StrainGroupingContext([context.project_path, context.database_path, context.fastq_path])
    sg = StrainGrouper(grouping_context)
    return sg.group()
            

def identify_sample_strains(strain_groups, context):
    print('\n\n=========== Identifying Strains ===========')
    write_group_header_to_report(context)
    identified_strains = []
    for group in strain_groups:
        identified_strains += characterise_group(group, context)

    return identified_strains


def write_group_header_to_report(context):
    output_filepath = context.project_name + '_detailed_report.tsv'
    with open(output_filepath, 'a') as fp:
        fp.write('strain name\tfilename\tgroup\tproportion within group\tMAPQ 60 count\tMAPQ 10 count\tMAPQ 2 count\n')


def characterise_group(group, context):
    print(f'\nGroup {group.id}')
    ga = GroupAnalyser(group, context)
    return ga.analyse_group()
    

def report_results(composition, context):
    composition = recalculate_sample_proportions(composition)
    print_results(composition)
    write_results(composition, context)


def recalculate_sample_proportions(composition):
    composition.sort(key=lambda x: x.sample_abundance, reverse=True)
    total_sample_characterised = sum([x.sample_abundance for x in composition])
    for strain in composition:
        strain.sample_abundance = strain.sample_abundance / total_sample_characterised * 100

    return composition


def print_results(composition):
    print('\n================================================== Final Composition ===================================================')
    print('{:>100}{:>20}'.format('strain name', 'sample abundance'))
    for strain in composition:
        print('{:>100}{:>20.2f}'.format(strain.name[:50], strain.sample_abundance))
    print()


def write_results(composition, context):
    output_filepath = context.project_name + '_brief_report.tsv'
    with open(output_filepath, 'w') as fp:
        fp.write('strain\tfilename\tsample proportion\n')
        for strain in composition:
            fp.write(f'{strain.name}\t{strain.filename}\t{round(strain.sample_abundance, 2)}\n')





if __name__ == '__main__':
    main(sys.argv[1:])


