import numpy as np
import math

def BoltsLoad(Fx,Fy,Fz,Mz,n,D2,e1,e2,e3,s2,t2,L):
    """ This function outputs an array with shear stresses
    for every bolt in the back plate"""
    # Author: Seppe
    rows = int(n/2)
    
    # Forces due to Fx
        # output in array
    Fxlist =[]
    for i in range(rows):
        Fxlist.append([((Fx*(s2+(t2/2)))/(2*((L/2)-e1)))/(n/2), -((Fx*(s2+(t2/2)))/(2*((L/2)-e1)))/(n/2)])
    Fxarray = np.array(Fxlist)

    # Forces due to Fy
        # output in array
    Fylist =[]
    for i in range(rows):
        Fylist.append([-Fy/n, -Fy/n])
    Fyarray = np.array(Fylist)

    # Forces due to Fz
        # Force closest to the neutral axis:
    denominatorfactor = 0
    Fzlist = []

    if (n%4) == 0:
        for i in range(int(n/4)):
            denominatorfactor += 2*i+1
        Fclosest = ((Fz*(s2+(t2/2)))/(2*e3*denominatorfactor))

        for i in range(int(n/2)-1,int(-(n/2)),-2):
            Fzlist.append([Fclosest*i, Fclosest*i])
    
    elif (n%4) != 0:
        for i in range(int((n-2)/4)):
            denominatorfactor += (i+1)^2
        Fclosest = ((Fz*(s2+(t2/2)))/(4*e3*denominatorfactor))

        for i in range(int((n-2)/4),int(-(((n-2)/4)+1)),-1):
            Fzlist.append([i*Fclosest, i*Fclosest])

    Fzarray = np.array(Fzlist)

    # Forces due to Mz
    FMzlist = []
    FMz1 = (-Mz)/(2*((L/2)-e1))
    FMz2 = (Mz)/(2*((L/2)-e1))

    for i in range(int(n/2)):
        FMzlist.append([FMz1,FMz2])

    FMzarray = np.array(FMzlist)

    # Sum of all forces per bolt
    Farray = Fxarray + Fyarray + Fzarray + FMzarray
        # Add all arrays
    
    # Surface Matrix
    Arealist = []
    Area = 2*math.pi*(D2/2)*t2
    for i in range(int((n/2))):
        Arealist.append([Area,Area])
    
    AreaArray = np.array(Arealist)

    # Shear Stress Matrix
    ShearStressArray = np.divide(Farray, AreaArray)

    return(ShearStressArray)


print(BoltsLoad(100,100,-100,100,10,10*10**(-3),30*10**(-3),30*10**(-3),30*10**(-3),50*10**(-3),5*10**(-3),200*10^(-3)))
