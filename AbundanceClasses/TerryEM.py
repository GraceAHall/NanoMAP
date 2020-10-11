
import json
from collections import defaultdict
import numpy as np
from copy import deepcopy




class TerryEM:
    def __init__(self, empirical_organism_summary, observed_counts):
        self.empirical_counts = empirical_organism_summary
        self.observed_counts = observed_counts
        self.create_objects()


    def create_objects(self):
        self.set_total_count()
        self.initialise_sample_proportions()
        self.set_basic_variables()
        self.create_transition_probabilities()
        self.create_count_matrix()


    def set_total_count(self):
        self.total_count = sum(self.observed_counts.values())


    def initialise_sample_proportions(self):
        num_organisms = len(self.empirical_counts)
        self.sample_proportions = defaultdict(int)
        for key in self.empirical_counts.keys():
            self.sample_proportions[key] = self.total_count / num_organisms


    def set_basic_variables(self):
        template = deepcopy(self.observed_counts)
        for key in template.keys():
            template[key] = 0
        self.predicted_counts = deepcopy(template)
        #self.residuals = deepcopy(template)
        self.colnames = list(self.observed_counts.keys())
        self.colnames.sort(key=lambda x: len(x))
        self.rownames = list(self.sample_proportions.keys())
        #self.max_ll = sum([self.observed_counts[col] * np.log(self.observed_counts[col]) for col in self.colnames])
        self.prev_ll = 0
        self.current_ll = 0
        self.likelihoods = []


    def create_transition_probabilities(self):
        self.transition_matrix = defaultdict(dict)
        for organism, counts in self.empirical_counts.items():
            total = sum(counts.values())
            for classification, count in counts.items():
                self.transition_matrix[organism][classification] = count / total


    def create_count_matrix(self):
        self.count_matrix = defaultdict(dict)
        for row in self.rownames:
            for col in self.colnames:
                self.count_matrix[row][col] = 0


    def iterate(self):
        self.iteration = 0
        self.update_predictions()
        #self.print_matrix(table='transition')
        #self.print_matrix(table='counts')
        self.update_parameters()

        while self.ll_delta > 0.1:
            self.iteration += 1
            print(f'\niteration {self.iteration}')
            print(f'log likelihood: {round(self.current_ll, 2)}')
            self.update_predictions()
            #self.print_matrix(table='counts')
            self.update_parameters()
                    
        self.set_final_proportions()
        self.print_final_proportions()
        #self.plot_ll()
        return self.sample_proportions, self.current_ll


    def update_predictions(self):
        self.update_count_matrix()
        self.predict_counts()
        

    def update_count_matrix(self):
        #prev_matrix = deepcopy(self.count_matrix)
        self.do_count_update()
        self.update_log_likelihood()
        self.update_ll_delta()
        #self.calculate_count_matrix_delta(prev_matrix)
        

    def do_count_update(self):
        sp = self.sample_proportions
        tm = self.transition_matrix
        cm = self.count_matrix

        for col in self.colnames:
            for row in self.rownames:
                cm[row][col] = sp[row] * tm[row][col]

    
    def update_log_likelihood(self):
        self.prev_ll = self.current_ll
        cm = self.count_matrix
        
        log_likelihood = 0
        for t in self.colnames:
            nt = self.observed_counts[t]
            qt = max(1, sum([cm[s][t] for s in self.rownames]))
            log_likelihood += nt * np.log(qt)
        
        self.current_ll = log_likelihood
        self.likelihoods.append([self.iteration, log_likelihood])


    def update_ll_delta(self):
        self.ll_delta = self.current_ll - self.prev_ll


    def calculate_count_matrix_delta(self, prev_matrix):
        total_delta = 0
        current_matrix = self.count_matrix

        for row in self.rownames:
            for col in self.colnames:
                total_delta += abs(current_matrix[row][col] - prev_matrix[row][col])
        self.delta = total_delta


    def predict_counts(self):
        cm = self.count_matrix
        for col in self.colnames:
            self.predicted_counts[col] = sum([cm[row][col] for row in self.rownames])


    def update_parameters(self):
        for s, m in self.sample_proportions.items():
            total = 0
            for t in self.colnames:
                nt = self.observed_counts[t]
                cst = self.count_matrix[s][t]
                pt = self.predicted_counts[t]
                if pt > 0:
                    total += nt * (cst / pt)
            self.sample_proportions[s] = total

    
    def print_matrix(self, table='norm_transition'):  
        if table == 'norm_transition':
            matrix = self.norm_transition_matrix
        elif table == 'transition':
            matrix = self.transition_matrix
        elif table == 'counts':
            matrix = self.count_matrix
        
        print('\n') 
        print(f'iteration: {self.iteration}')    
        print('{:^10}{:^10}'.format('organism', 'Ps'), end='')
        for col in self.colnames:
            print('{:^10}'.format(col), end='')
        print()

        for row in self.rownames:
            current_proportion = self.sample_proportions[row]
            matrix_row = matrix[row]
            if matrix == self.count_matrix:
                print('{:^10}{:^10}'.format(row, int(current_proportion)), end='')
                for col in self.colnames:
                    print('{:^10}'.format(int(matrix_row[col])), end='')
                print()
            else:
                print('{:^10}{:^10.3f}'.format(row, int(current_proportion)), end='')
                for col in self.colnames:
                    print('{:^10.3f}'.format(matrix_row[col]), end='')
                print()
        
        print('{:^10}{:^10}'.format('', 'predicted'), end='')
        for col in self.colnames:
            print('{:^10}'.format(int(self.predicted_counts[col])), end='')
        print()
        
        print('{:^10}{:^10}'.format('', 'observed'), end='')
        for col in self.colnames:
            print('{:^10}'.format(self.observed_counts[col]), end='')
        print()

        print('previous log likelihood: {:0.3f}'.format(self.prev_ll))
        print('current log likelihood: {:0.3f}'.format(self.current_ll))
        print('delta: {:0.3f}'.format(self.ll_delta))


    def set_final_proportions(self):
        for symbol, count in self.sample_proportions.items():
            self.sample_proportions[symbol] = count / self.total_count


    def print_final_proportions(self):
        print('\n\n### sample proportions ###')
        print('{:^15}{:^15}'.format('organism', 'proportion'))
        for symbol, count in self.sample_proportions.items():
            print('{:^15}{:^15.2f}'.format(symbol, count))

