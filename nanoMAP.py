

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
    check_environment()
    opts, args = getopt.getopt(argv, "ht:r:d:p:", ["map-pb", "map-ont", "no-initial-alignment", "allow-extrachromosomal"])
    context = ProgramContext(opts) 

    params = [opt[0] for opt in opts]
    if "--no-initial-alignment" in params:
        perform_noalign_setup(context)
    else:
        perform_full_setup(context)
        initial_alignment(context)

    strain_groups = group_strains(context)
    identified_strains = identify_sample_strains(strain_groups, context)
    composition = composition_analysis(identified_strains, context)
    report_results(composition, context)
                           

def check_environment():
    if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 6):
        raise Exception("Must be using Python 3.6 or greater")


def perform_full_setup(context):
    fs = FileSetup(context)
    fs.full_setup()


def perform_noalign_setup(context):
    fs = FileSetup(context)
    fs.noalign_setup()


def initial_alignment(context):
    mm2 = Minimap2Runner(context)
    mm2.initial_alignment()
    

def group_strains(context):
    grouping_context = StrainGroupingContext([context.project_path, context.database_path, context.fastq_path])
    sg = StrainGrouper(grouping_context)
    return sg.group()
            

def identify_sample_strains(strain_groups, context):
    identified_strains = []
    for group in strain_groups:
        identified_strains += characterise_group(group, context)

    return identified_strains


def characterise_group(group, context):
    ga = GroupAnalyser(group, context)
    return ga.analyse_group()
    

def composition_analysis(identified_strains, context):
    cpa = CompositionAnalyser(identified_strains, context)
    return cpa.analyse()
 

def report_results(composition, context):
    composition.sort(key=lambda x: x.sample_abundance, reverse=True)
    print('final composition:')
    with open(context.project_path + '/report.tsv', 'w') as fp:
        for strain in composition:
            print('{:>50}{:>15.2f}'.format(strain.name[:50], strain.sample_abundance))
            fp.write(f'{strain.name}\t{",".join(strain.accessions)}\t{round(strain.sample_abundance, 2)}\n')


if __name__ == '__main__':
    main(sys.argv[1:])


