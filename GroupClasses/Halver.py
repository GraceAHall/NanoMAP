

import math


class Halver:
    def __init__(self, characterisation, context):
        self.characterisation = characterisation
        self.min_abundance = 0.05 * (100 / characterisation[0].group_abundance) 
        #print(self.min_abundance)


    def halve(self):
        # top 75 % mapq1 
        # top 66 % collinearity
        characterisation = self.characterisation

        #print(self.min_abundance)

        strain_quota = math.ceil(0.5 * len(characterisation))
        remaining_strains = math.ceil(0.5 * len(characterisation))
        characterisation = [strain for strain in characterisation if strain.naive_abundance > self.min_abundance]

        if len(characterisation) < strain_quota:
            return characterisation

        # iteratively place strains in new_characterisation
        new_characterisation = []
        
        print(f'\nselecting {strain_quota} strains\n')
        # put top things 
        characterisation.sort(key=lambda x: x.mapq_dict[60], reverse=True)
        for i in range(min(2, remaining_strains)):
            strains_to_delete = []

            if characterisation[i].mapq_dict[60] > 2:
                new_characterisation.append(characterisation[i])
                strains_to_delete.append(i)

            for i in reversed(strains_to_delete):
                del characterisation[i]

        remaining_strains = strain_quota - len(new_characterisation)

        if len(new_characterisation) >= strain_quota:
            return new_characterisation

        characterisation.sort(key=lambda x: x.mapq_dict[10], reverse=True)
        for i in range(min(2, remaining_strains)):
            strains_to_delete = []
            if characterisation[i].mapq_dict[10] > 2:
                new_characterisation.append(characterisation[i])
                strains_to_delete.append(i)
            for i in reversed(strains_to_delete):
                del characterisation[i]

        remaining_strains = strain_quota - len(new_characterisation)
        
        if len(new_characterisation) >= strain_quota:
            return new_characterisation

        characterisation.sort(key=lambda x: x.mapq_dict[2], reverse=True)
        new_characterisation += characterisation[:remaining_strains]

        return new_characterisation
    