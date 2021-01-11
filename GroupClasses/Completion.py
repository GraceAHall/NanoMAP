

 




class CompletionJudge:
    def __init__(self, characterisation): 
        self.characterisation = characterisation
        self.min_num_strains = 4
        self.min_cumulative_mapq_60 = 1000


    def judge(self):
        characterisation = self.characterisation
        if len(characterisation) <= self.min_num_strains:
            return True

        characterisation.sort(key=lambda x: x.mapq_dict[60], reverse=True)
        if characterisation[0].mapq_dict[60] > 100 and characterisation[1].mapq_dict[60] < 10:
            return True

        return False