


class GroupPruner:
    def __init__(self, groups):
        self.groups = groups
        self.min_group_bases = 1 * 1000000  # 1mb sequence data


    def prune(self):
        self.format_groups_to_list()
        #self.remove_single_strain_groups()
        self.remove_low_abundance_groups()
        return self.groups


    def format_groups_to_list(self):
        self.groups = list(self.groups.values())
        self.groups.sort(key=lambda x: x.abundance, reverse=True)


    def remove_single_strain_groups(self):
        groups_to_delete = []

        for i, group in enumerate(self.groups):
            if len(group.identifiers) < 2:
                groups_to_delete.append(i)
        
        for index in reversed(groups_to_delete):
            del self.groups[index]


    def remove_low_abundance_groups(self):
        groups_to_delete = []
        
        for i, group in enumerate(self.groups):
            if group.basecount < self.min_group_bases:
                groups_to_delete.append(i)
        
        for index in reversed(groups_to_delete):
            del self.groups[index]
