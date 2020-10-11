




class StrainPicker:
    def __init__(self, characterisation):
        self.characterisation = characterisation


    def pick(self):
        if len(self.characterisation) <= 1:
            return self.characterisation

        characterisation = self.characterisation
        characterisation.sort(key=lambda x: x.mapq_dict[60], reverse=True)

        # just avoiding 50 * 0
        tophit = characterisation[0].mapq_dict[60]
        secondhit = max(1, characterisation[1].mapq_dict[60])

        # checking the 10x condition
        if tophit >= 10 * secondhit:
            characterisation = [characterisation[0]]

        else:
            characterisation = [strain for strain in characterisation if strain.mapq_dict[60] > 0]
            characterisation.sort(key=lambda x: x.mapq_dict[10], reverse=True)
            characterisation = characterisation[:5]

        return characterisation


