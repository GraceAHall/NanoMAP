
from modules.io import load_json, write_reads_from_dict, write_group_genomes
from GeneralClasses.Strain import Strain, Contig 
from GroupClasses.Halver import Halver


from collections import defaultdict


class Characteriser:
    def __init__(self, classifications, collinearities, group_id, group_abundance, context):
        self.group_id = group_id
        self.group_abundance = group_abundance
        self.classifications = classifications
        self.collinearities = collinearities
        self.context = context
        self.min_cec_ratio = 2

        self.af_dict = load_json(context.database_path + 'taxonomy/accessions_filenames.json')
        self.ag_dict = load_json(context.database_path + 'taxonomy/accessions_genomes.json')


    def format_characterisation(self):
        # initialise characterisation
        characterisation = self.initialise_characterisation()
        characterisation = self.add_basecounts(characterisation)
        characterisation = self.add_collinearities(characterisation)
        characterisation = self.remove_extrachromosomal_dominant_strains(characterisation)
        
        characterisation = self.add_mapqs(characterisation)
        characterisation = self.add_abundances(characterisation)
        characterisation = self.flatten_characterisation(characterisation)
        characterisation = self.set_primary_names(characterisation)
        return characterisation


    def initialise_characterisation(self):
        characterisation = {}    
        classifications = self.classifications    
        af_dict = self.af_dict
        ag_dict = self.ag_dict

        accessions = []
        for alignments in classifications.values():
            accessions += [al[1] for al in alignments]
        accessions = set(accessions)

        for accession in accessions:
            accession = accession.upper()
            try:
                filename = af_dict[accession]
                contig_name = ag_dict[accession]
            except KeyError:
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
        classifications = self.classifications
        af_dict = self.af_dict
        ag_dict = self.ag_dict

        for read_id, alignments in classifications.items():
            equal_bases_share = alignments[0][5] / len(alignments)
            for al in alignments:
                accession = al[1].upper()
                try:
                    filename = af_dict[accession]
                    contig_name = ag_dict[accession]
                except KeyError:
                    print(f'keyerror at add_abundances: {accession}')
                    filename = accession
                    contig_name = accession

                characterisation[filename].contigs[accession].basecount += equal_bases_share

        for strain in characterisation.values():
            strain.set_basecount()

        return characterisation


    def add_collinearities(self, characterisation):
        collinearities = self.collinearities
        taxonomy = self.af_dict

        for filename, strain in characterisation.items():
            strain.collinearity = collinearities[filename]

        return characterisation


    def add_mapqs(self, characterisation):
        classifications = self.classifications
        ag_dict = self.ag_dict
        af_dict = self.af_dict
        
        for read_id, alignments in classifications.items():
            for al in alignments:
                try:
                    filename = af_dict[al[1].upper()]
                    characterisation[filename].mapq_alignments.append(al[7])
                except KeyError:
                    if al[1].upper() in ag_dict:
                        ident = ag_dict[al[1].upper()]
                    else:
                        ident = al[1]
                    print('taxonomy error. is the database build current?')
                    continue 

        for strain in characterisation.values():
            if len(strain.mapq_alignments) > 0:
                mapq_dict = defaultdict(int)
                for mapq_score in strain.mapq_alignments:
                    if mapq_score >= 1:
                        mapq_dict[1] += 1
                    if mapq_score >= 2:
                        mapq_dict[2] += 1
                    if mapq_score >= 5:
                        mapq_dict[5] += 1
                    if mapq_score >= 10:
                        mapq_dict[10] += 1
                    if mapq_score >= 60:
                        mapq_dict[60] += 1
                strain.mapq_dict = mapq_dict

        return characterisation


    def remove_extrachromosomal_dominant_strains(self, characterisation):
        strains_to_delete = []

        for filename, strain in characterisation.items():
            strain.set_chromosomal_extrachromosomal_ratio()
            if strain.cec_ratio < self.min_cec_ratio:
                strains_to_delete.append(filename)
        
        print(f'{len(strains_to_delete)} strains are extrachromosomal dominant')
        for filename in reversed(strains_to_delete):
            del characterisation[filename]
        
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



