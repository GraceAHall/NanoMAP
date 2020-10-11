

from CompositionClasses.PafProcessing import PafProcessor
from CompositionClasses.Characterisation import Characteriser
from AbundanceClasses.ProportionEstimator import ProportionEstimator


import subprocess
import os


class CompositionAnalyser:
    def __init__(self, characterisation, context):
        self.characterisation = characterisation
        self.context = context


    def analyse(self):
        self.estimate_abundances()
        return self.identified_strains


    def estimate_abundances(self):
        pe = ProportionEstimator(self.characterisation, self.context)
        self.identified_strains = pe.estimate()



