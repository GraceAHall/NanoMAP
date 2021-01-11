

import math


class Halver:
    def __init__(self, characterisation, context):
        self.old_characterisation = characterisation
        self.new_characterisation = []
        self.strains_to_delete = []

        # amount DNA attributed to strain needed for strain consideration
        self.min_abundance = 0.05 * (100 / characterisation[0].group_abundance) 

        # number of strains to keep after selection process
        self.quota = math.ceil(0.5 * len(characterisation))
        


    def halve(self):
        print(f'\nselecting {self.quota} strains\n')
        self.remove_low_abundance_strains()

        for score in [60, 10, 2]:
            if self.has_reached_quota():
                return self.new_characterisation
            self.add_strains_by_mapq_score(score)
            self.delete_strains_from_old_characterisation()
            
        return self.new_characterisation
            

    def remove_low_abundance_strains(self):
        self.old_characterisation = [s for s in self.old_characterisation if s.naive_abundance > self.min_abundance]


    def has_reached_quota(self):
        if len(self.new_characterisation) >= self.quota:
            return True
        return False


    def add_strains_by_mapq_score(self, score):
        self.old_characterisation.sort(key=lambda x: x.mapq_dict[60], reverse=True)

        i = 0
        while len(self.new_characterisation) < self.quota and i < len(self.old_characterisation):
            strain = self.old_characterisation[i]
            if strain.mapq_dict[score] > 2:
                self.new_characterisation.append(strain)
                self.strains_to_delete.append(i)
            i += 1


    def delete_strains_from_old_characterisation(self):
        for i in reversed(self.strains_to_delete):
            del self.old_characterisation[i] 
        self.strains_to_delete = []


    
  
    
    

