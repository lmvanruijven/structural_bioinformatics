#!/usr/bin/env python
#HHPred

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='alignment_PSA.ali',
              knowns='3zih', sequence='T0868',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
