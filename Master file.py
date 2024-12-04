#This is the master file 
#Please do not edit this unless you're doing something useful
import numpy as np
import math
import matplotlib.pyplot as plt
import FailureModeFunctions

#Import files
import constants as cst
import Configuration_I as config

# Perform checks
    # => Safety factors

# Backplate design

for D2 in range(5):
    for L in range(5):
        for t2 in range(5):
            for n in range(5):
                FailureModeFunctions.PullThrough()