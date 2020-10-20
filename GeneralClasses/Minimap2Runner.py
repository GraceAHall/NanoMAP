

import os
import subprocess



class Minimap2Runner:
    def __init__(self, context):
        self.context = context

    
    def initial_alignment(self):
        print('running minimap2 on full database')
        ct = self.context

        paf_output = ct.project_path + '/runtimefiles/pafs/full_alignment.paf'
        preset = ct.read_technology
        threads = ct.threads
        index = ct.database_path + 'database.mmi'
        reads = ct.fastq_path
        
        with open(paf_output, 'w') as outfile:
            subprocess.run(['minimap2', '-x', preset, '-I', str(ct.max_memory) + 'G', '-t', threads, index, reads], stdout=outfile)



