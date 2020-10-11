



class GroupAbundanceCalculator:
    def __init__(self, strain_groups, group_membership, strain_basecounts):
        self.strain_groups = strain_groups
        self.group_membership = group_membership
        self.strain_basecounts = strain_basecounts


    def calculate(self):
        self.set_group_base_counts()
        self.set_strain_base_counts()
        self.set_group_abundances()
        self.print_group_abundances()
        return self.strain_groups


    def set_group_base_counts(self):
        strain_basecounts = self.strain_basecounts
        group_membership = self.group_membership

        for identifier, basecount in strain_basecounts.items():
            try:
                group_id = group_membership[identifier]
                self.strain_groups[group_id].basecount += basecount
            except KeyError:
                pass


    def set_strain_base_counts(self):
        strain_basecounts = self.strain_basecounts
        group_membership = self.group_membership

        for identifier, basecount in strain_basecounts.items():
            try:
                group_id = group_membership[identifier]
                self.strain_groups[group_id].strain_basecounts[identifier] = basecount
            except KeyError:
                pass
        

    def set_group_abundances(self):
        total_basecount = sum([x.basecount for x in self.strain_groups.values()])
        for group in self.strain_groups.values():
            group.abundance = group.basecount / total_basecount * 100


    def print_group_abundances(self):
        print('\ngroup abundances:')
        for group in list(self.strain_groups.values())[:20]:
            print('{:10}{:>10.1f}'.format(group.id, group.abundance))