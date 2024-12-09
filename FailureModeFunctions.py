import math
import numpy as np
import Configuration_I as c
import csv
#

def Forces(Fx, Fy, Fz, H, W):
    # Author: Seppe
    """
    Returns a list [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz]
    Fx = Force in x-direction
    Fy = Force in y-direction
    Fz = Force in z-direction
    H  = Height of the solar panel in z-direction
    W  = length of the solar panel in y-direction
    """
    #Forces
    Az = -Fz/2 #why negative
    Bz = -Fz/2
    Ax = -Fx/2
    Bx = -Fx/2
    Ay = (Fz*W-Fy*H)/(2*H)
    By = -Fy-Ay


    # Moments
    M_Ay = -(Fx * H)/8
    M_By = (Fx * H)/8
    M_Az = (Fx * W)/4
    M_Bz = (Fx * W)/4
    
    
    My = max(abs(M_Ay), abs(M_By))  
    Mz = max(abs(M_Az), abs(M_Bz))
    Fortaz = abs(Az)
    Fortax = abs(Ax)

    return [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz, My, Mz, Fortax, Fortaz]

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

def FlangeFailure(W,D,t,S_ty,F_y,F_z):
    A_br = D*t
    A_t = (W-D)*t

    # K_bry, fig. D1.14 in Bruhn, taken max curve with out of bounds if the relevant curve diverges
    if (0.5*W)/D < 3.3*(t/D)+0.57:
       K_bry = -1.25+2.98*((0.5*W)/D)-0.998*((0.5*W)/D)**(2)+0.11*((0.5*W)/D)**3
    else:
        # print("e/D conbined with D/t is out of bounds")
        K_bry = 100000000

    # K_ty, fig. D1.15 in Bruhn, only curve 3
    A_2 = 0.5*(W-D)*t
    A_1 = (A_2+0.5*D-0.5*D*math.cos(math.radians(45)))*t
    A_av = 6/((4/A_1)+(2/A_2))
    K_ty = -4.72*10**(-3)+1.39*(A_av/A_br)-0.341*(A_av/A_br)**2

    # K_t, fig. D1.12 in Bruhn, only linear part otherwise out of bounds
    if W/D <= 2.9: 
        K_t = (W/D)*(-0.08/1.4)+1.0857
    else:
        K_t = 100000000
        # print("W/D out of bounds (over 2.9)")

    # transverse
    P_ty = K_ty*A_br*S_ty
    R_tr = abs(F_z)/P_ty

    # axial
    P_y = K_t*S_ty*A_t
    P_bry = K_bry*S_ty*A_br

    if P_y > P_bry:
        R_a = abs(F_y)/P_bry
    else:
        R_a = abs(F_y)/P_y

    SF = (1/((R_a**1.6+R_tr**1.6)**0.625))-1
    print("K_bry: ", K_bry, " K_ty: ", K_ty, " K_t: ", K_t, " P_ty: ", P_ty, " P_y: ", P_y, " P_bry: ", P_bry, " F_y: ", F_y, " F_z: ", F_z)
    return SF




def BearingFailure(Ax, Az, My, D_2, t_2, L, n):
    x = L/2 - 1.5*D_2

    if (n == 4):
        z = 1.25 * D_2
        F_xbolt = Ax / n + My / (n * x)
        F_zbolt = Az / n + My / (n * z)
    if (n == 6):
        z = 2.5 * D_2
        F_xbolt = Ax / n + My / (n * x)
        F_zbolt = Az / n + My / ((n - 2) * z)
        
    if (n == 8):
        z = 3.75 * D_2
        F_xbolt = Ax / n + My / (n * x)
        F_zbolt = Az / n + My / ((n - 4) * z)

    P = (F_xbolt ** 2 + F_zbolt ** 2) ** 0.5

    sigma = 1.2 * P / (D_2 * t_2)  # max stress experienced by the bolt
    # then compare sigma to the one of the material max strenght and see how to lighten up the hinge

    sigmamaterial = c.bolt.shear_strength   

    return sigma / sigmamaterial


def Thermal(t_1, materialid, D_2, c.materials, c.bolt)
    sa = 4* t_1 /(c.materials[materialid].elastic_modulus*3.14*1.5*D_2**2)
    sb = 1/ c.bolt.elastic_modulus * 4* t_1 /(3.14* D_2**2/4)
    
    stupidletter = sa/(sa+sb)
    
    stress = (c.materials[materialid].thermal_coef- c.bolt.thermal_coef)*250 * c.bolt.elastic_modulus *3.14 * D_2**2/4 * (1-stupidletter)
    
    return stress
    
    



def MassCalc(s2, D1, t1, w, t2, L, n, D2, rho):
    VolumeBP = t2*w*L - n*D2*t2
    VolumeF  = s2*w*t1 + 0.5*math.pi*((0.5*w)**2)*t1 - math.pi*((D1/2)**2)*t1

    Volume   = VolumeBP + 2 * VolumeF

    mass = Volume * rho

    return mass

# fortele=Forces(c.Fx, c.Fy, c.Fz, c.H, c.W)
# Fortax=fortele[12]
# Fortaz=fortele[13]
# My=fortele[10]

# LOOP FOR LUG A:
# Author: Seppe
# Enter general forces here: (Taken from the report)
Fx = 215.82
Fy = 215.82
Fz = -686.7

# Then the forces at lug A are calculated:
Ax = Forces(Fx,Fy,Fz,0.450,0.975)[0]
Ay = Forces(Fx,Fy,Fz,0.450,0.975)[1]
Az = Forces(Fx,Fy,Fz,0.450,0.975)[2]

M_Ay = Forces(Fx,Fy,Fz,0.450,0.975)[6]
M_Az = Forces(Fx,Fy,Fz,0.450,0.975)[7]

# Then we also need some material properties:
Tau_max = 27 *10**(9) # Pa
S_ty = 10
rho = 2000

# First iterate over flanges to determine s2

data = [['Iteration', 'D1','D2', 'L', 'W', 't1', 't2', 'n', 'SF Pullthrough', 'SF Flangefailure', 'SF Bearing', 'mass']]

iteration = 1

<<<<<<< HEAD
for D2 in np.arange(0.001,0.02,0.001):
    for L in np.arange(0.01,0.2,0.01):
=======
for D2 in np.arange(0.002,0.010,0.002):
    for L in np.arange(0.02,0.3,0.02):
>>>>>>> 4bd0384f4aed0a855b5eb41df4db99c3dd4a85c4
        for t2 in np.arange(0.001,0.01,0.001):
            for n in range(4,10,2):
<<<<<<< HEAD
                for D1 in np.arange(0.060,0.080,2):
=======
                e1 = 1.5 * D2
                e3 = 2.5 * D2
                W = 2* e1 + (n/2)*D2 + ((n/2)-1)* e3
                for D1 in np.arange(0.005,W-0.001,0.005):
>>>>>>> 4bd0384f4aed0a855b5eb41df4db99c3dd4a85c4
                    Pullthrougharray = Pullthrough(Ax, Ay, Az, M_Az, n, D2, 1.5*D2, 1.5*D2, 2.5*D2, 2*D1, t2, L)
                    Tau_max_list = []
                    for i in range(int((n/2))):
                        Tau_max_list.append([Tau_max, Tau_max])
                    Tau_max_array = np.array(Tau_max_list)
                    AbsPullthrougharray = abs(Pullthrougharray)
                    Marginsarray = Tau_max_array - AbsPullthrougharray
                    negative = np.any(Marginsarray<0)
                    if negative == False:
                        max = np.max(AbsPullthrougharray)
                        SFPullthrough = Tau_max / max
                        

                    e1 = 1.5 * D2
                    e3 = 2.5 * D2

                    W = 2* e1 + (n/2)*D2 + ((n/2)-1)* e3

                
<<<<<<< HEAD
                    for t1 in np.arange(0.001,0.011,0.001):
=======
                    for t1 in np.arange(0.003,0.015,0.003):
>>>>>>> 4bd0384f4aed0a855b5eb41df4db99c3dd4a85c4
                        SFflange = FlangeFailure(W,D1,t1,S_ty, Ay, Az)
                        SFbearing = BearingFailure(Ax,Az,M_Ay,D2,t2,n, )

                        mass = MassCalc(2*D1, D1, t1, W, t2, L,n,D2, rho)

                        data.append([iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass])
                        print("Iteration: ", iteration)
                        iteration += 1

with open('Designpoints.csv', 'w', newline = '') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
    csvfile.close()