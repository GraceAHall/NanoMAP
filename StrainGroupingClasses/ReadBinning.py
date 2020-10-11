



class ReadBinner:
    def __init__(self, groups, group_membership, read_classifications):
        self.groups = groups
        self.group_membership = group_membership
        self.read_classifications = read_classifications
    
    
    def bin(self):
        groups = self.groups
        group_membership = self.group_membership
        read_classifications = self.read_classifications

        binned = 0
        unbinned = 0
        for read in read_classifications:
            group_assignments = set()
            classifications = read.classifications
            for identifier in classifications:
                group_id = group_membership[identifier]
                group_assignments.add(group_id)
            
            if len(group_assignments) == 1:
                binned += 1
                group_id = group_assignments.pop()
                groups[group_id].read_ids.add(read.read_id)
            else:
                unbinned += 1

        return groups