
from modules.io import load_json
from GeneralClasses.Characterisation import SampleCharacterisation



class SampleCharacteriser:
    def __init__(self, reads_alignments, context):
        self.reads_alignments = reads_alignments
        self.context = context

        self.af_dict = load_json(self.context.database_path + 'taxonomy/accessions_filenames.json')
        self.an_dict = load_json(self.context.database_path + 'taxonomy/accessions_genomes.json')


    def characterise(self):
        self.fill_characterisation_structure()
        self.characterisation.reformat_strains_to_lists()
        self.characterisation.set_sampled_base_count()
        self.characterisation.set_sample_abundances()
        self.characterisation.remove_extrachromosomal_dominant_strains()
        return self.characterisation


    def fill_characterisation_structure(self):
        sc = SampleCharacterisation(self.af_dict, self.an_dict)
        reads_alignments = self.reads_alignments
        for read_id, alignments in reads_alignments.items():
            sc.update_characterisation(read_id, alignments)
        self.characterisation = sc
        # TODO I AM HERE 
     
        


        
