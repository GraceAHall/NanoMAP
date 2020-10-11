


from modules.io import load_json, write_reads_from_dict, write_group_genomes
from GeneralClasses.Strain import Strain, Contig 

#from GroupClasses.PafProcessor import PafProcessorGroup
#from GroupClasses.StrainPicker import StrainPicker
#from GeneralClasses.Characterisation import SampleCharacterisation

from collections import defaultdict


class Characteriser:
    def __init__(self, best_alignments, context):
        self.best_alignments = best_alignments
        self.context = context

        self.af_dict = load_json(context.database_path + 'taxonomy/accessions_filenames.json')
        self.ag_dict = load_json(context.database_path + 'taxonomy/accessions_genomes.json')


    def characterise(self):
        characterisation = self.initialise_characterisation()
        characterisation = self.add_basecounts(characterisation)
        characterisation = self.add_abundances(characterisation)
        characterisation = self.flatten_characterisation(characterisation)
        characterisation = self.set_primary_names(characterisation)
        return characterisation


    def initialise_characterisation(self):
        characterisation = {}    
        alignments = self.best_alignments    
        af_dict = self.af_dict
        ag_dict = self.ag_dict

        accessions = [al[1] for al in alignments]

        for accession in accessions:
            accession = accession.upper()
            try:
                filename = af_dict[accession]
                contig_name = ag_dict[accession]
            except KeyError:
                print(f'cant find filename/contigname for {accession}')
                filename = accession
                contig_name = accession

            if filename not in characterisation:
                new_strain = Strain(filename)
                new_strain.group_id = self.group_id
                new_strain.group_abundance = self.group_abundance
                characterisation[filename] = new_strain

            characterisation[filename].contigs[accession] = Contig(accession, contig_name)

        return characterisation


    def add_basecounts(self, characterisation):
        alignments = self.best_alignments
        af_dict = self.af_dict
        ag_dict = self.ag_dict

        for al in alignments:
            accession = al[1].upper()
            read_length = al[5]
            try:
                filename = af_dict[accession]
                contig_name = ag_dict[accession]
            except KeyError:
                print(f'keyerror at add_abundances: {accession}')
                filename = accession
                contig_name = accession

            characterisation[filename].contigs[accession].basecount += read_length

        for strain in characterisation.values():
            strain.set_basecount()

        return characterisation


    def add_abundances(self, characterisation):
        total_basecount = sum([strain.basecount for strain in characterisation.values()])
        for strain in characterisation.values():
            strain.naive_abundance = strain.basecount / total_basecount * 100

        return characterisation


    def flatten_characterisation(self, characterisation):
        characterisation = list(characterisation.values())
        characterisation.sort(key=lambda x: x.naive_abundance, reverse=True)
        for strain in characterisation:
            strain.accessions = list(strain.contigs.keys())
        return characterisation


    def set_primary_names(self, characterisation):
        ag_dict = self.ag_dict

        for strain in characterisation:
            for accession in strain.contigs:
                contig_name = ag_dict[accession.upper()]
                if 'plasmid' not in contig_name and 'mitochondr' not in contig_name:
                    contig_name = contig_name.rstrip(', complete genome')
                    contig_name = contig_name.rstrip(', complete sequence')
                    contig_name = contig_name.rstrip(' chromosome')
                    strain.name = contig_name
                    break
        
        return characterisation





