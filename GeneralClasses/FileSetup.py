



import os
import shutil
import sys
import subprocess

class FileSetup:
    def __init__(self, context):
        self.context = context


    def full_setup(self):
        self.check_if_project_exists()
        self.make_project_folder()
        self.make_project_subfolders()


    def check_if_project_exists(self):
        pp = self.context.project_path
        if os.path.exists(pp):
            project_name = pp.lstrip('projects/')
            confirmation = input(f'the project {project_name} already exists - would you like to wipe it and continue? [y/n]')
            if confirmation == 'y' or confirmation == 'Y':
                shutil.rmtree(pp, ignore_errors=True)
            else:
                print('exiting')
                sys.exit()


    def make_project_folder(self): 
        pp = self.context.project_path       
        try:
            os.mkdir('projects')
        except FileExistsError:
            pass
    
        os.mkdir(pp)


    def make_project_subfolders(self):
        pp = self.context.project_path
        os.mkdir(pp + '/runtimefiles')
        os.mkdir(pp + '/runtimefiles/pafs')
        os.mkdir(pp + '/runtimefiles/group_databases')
        os.mkdir(pp + '/runtimefiles/group_fastqs')
        os.mkdir(pp + '/runtimefiles/characterisations')
        os.mkdir(pp + '/runtimefiles/abudance_estimation')
        open(pp + '/banlist.txt', 'a') 


    def noalign_setup(self):
        pp = self.context.project_path
        subprocess.call(f'rm {pp}/runtimefiles/group_databases/*', shell=True)
        subprocess.call(f'rm {pp}/runtimefiles/group_fastqs/*', shell=True)
        # more here please

    

