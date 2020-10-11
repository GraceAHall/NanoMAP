

from StrainGroupingClasses.StrainGroup import StrainGroup



class Grouper:
    def __init__(self, read_counts, co_occurance):
        self.co_occurance = co_occurance
        self.read_counts = list(read_counts.items())
        self.read_counts.sort(key=lambda x: x[1], reverse=True)
        self.min_co_occurance_rate = 0.5


    def group(self):
        self.form_groups()
        self.set_group_membership()
        return self.groups, self.group_membership


    def form_groups(self):
        read_counts = self.read_counts
        co_occurance = self.co_occurance

        min_rate = self.min_co_occurance_rate

        strain_groups = {}
        group_id = 0
        while len(read_counts) > 0:
            identifier1, _ = read_counts.pop(0)
            group = StrainGroup(group_id)
            group.identifiers.add(identifier1)
            items_to_delete = []

            i = 0
            for identifier2, __ in read_counts:
                rate = co_occurance[identifier2][identifier1]
                if rate > min_rate:
                    group.identifiers.add(identifier2)
                    items_to_delete.append(i)
                i += 1
            
            strain_groups[group.id] = group

            for index in reversed(items_to_delete):
                del read_counts[index]
            
            group_id += 1
        
        self.groups = strain_groups

    
    def set_group_membership(self):
        group_membership = {}
        for group_id, group in self.groups.items():
            for identifier in group.identifiers:
                group_membership[identifier] = group.id
        self.group_membership = group_membership








