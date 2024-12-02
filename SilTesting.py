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





# Sil testing
SF = FM.FlangeFailure(config.w, config.D1, config.t1, 70000000, Forces[2], Forces[3])
D_range = 4 * 10**(-3)
W_range = 4 * 10**(-3)
t_range = 1 * 10**(-3)

D_step = 0.25 * 10**(-3)
W_step = 0.25 * 10**(-3)
t_step = 0.1 * 10**(-3)

LSF = 100

Dx = config.D1-D_range
while Dx < (config.D1+D_range):
    Wx = config.w-W_range
    while Wx < (config.w+W_range):
        tx = config.t1-t_range
        while tx < (config.w+W_range):
            SF = FM.FlangeFailure(Wx, Dx, tx, 70000000, Forces[2], Forces[3])
            print(Dx, Wx, tx, SF)
            if SF < LSF and SF > 1.25:
                LSF = SF
                Config = [Dx, Wx, tx]
            tx = tx + t_step
        Wx = Wx + W_step
    Dx = Dx + D_step
print(Config)

print("SF: ", FM.FlangeFailure(Config[1], Config[0], Config[2], 70000000, Forces[2], Forces[3]))

