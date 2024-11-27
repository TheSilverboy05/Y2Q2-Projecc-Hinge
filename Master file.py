#This is the master file 
#Please do not edit this unless you're doing something useful
import numpy as np
import math
import matplotlib.pyplot as plt
import FailureModeFunctions as FM

#Import files
import constants as cst
import Configuration_I as config

# Perform checks

Forces = FM.Forces(config.Fx, config.Fy, config.Fz, config.H, config.W)

SF1 = FM.FlangeFailure()