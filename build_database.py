

import sys
import getopt
import os
from GeneralClasses.Context import DatabaseBuildingContext
from DatabaseClasses.IndexBuilder import IndexBuilder
from DatabaseClasses.TaxonomyBuilder import TaxonomyBuilder



# takes the folder to process as positional arg 
def main(argv): 
    opts, args = getopt.getopt(argv, "hd:t:p:m:", ["taxonomy-only", "database-only", "rebuild"])
    context = DatabaseBuildingContext(opts)

    check_existing_build(context)

    for opt, arg in opts:
        if opt == '--taxonomy-only':
            build_taxonomy(context)
            sys.exit()
        elif opt == '--database-only':
            build_minimap2_index(context)
            sys.exit()
        
    build_taxonomy(context)
    build_minimap2_index(context)
    print('\nfinished database build.')


def check_existing_build(context):
    path = context.database_path
    if os.path.exists(path + 'taxonomy/accessions_genomes.json') or os.path.exists(path + 'database.mmi'):
        if context.rebuild is False:
            print('database has already been built. If you wish to override, give the "--rebuild" option to this program')
            sys.exit()


def build_minimap2_index(context):
    ib = IndexBuilder(context)
    ib.build()


def build_taxonomy(context):
    tb = TaxonomyBuilder(context)
    tb.build()


if __name__ == '__main__':
    main(sys.argv[1:])








    






