import math
import numpy as np
import Configuration_I as c
import csv
#

def Forces(Fx, Fy, Fz, H, S):
    # Author: Seppe
    """
    Returns a list [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz]
    Fx = Force in x-direction
    Fy = Force in y-direction
    Fz = Force in z-direction
    H  = Height of the solar panel in z-direction
    S  = length of the solar panel in y-direction
    """
    #Forces
    Az = -Fz/2 #why negative
    Bz = -Fz/2
    Ax = -Fx/2
    Bx = -Fx/2
    Ay = -Fy/2+Fz*S/H
    By = -Fy-Ay


    # Moments
    M_Ay = -(Fx * H)/8
    M_By = (Fx * H)/8
    M_Az = (Fx * S)/2
    M_Bz = (Fx * S)/2
    
    
    My = max(abs(M_Ay), abs(M_By))  
    Mz = max(abs(M_Az), abs(M_Bz))
    Fortaz = abs(Az)
    Fortax = abs(Ax)

    return [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz, My, Mz, Fortax, Fortaz]

def FlangeLoads(Ax,Ay,Az,M_Ay,M_Az,h,t1):
    Rx=-Ax
    Ry=-Ay
    Rz=-Az
    My=-M_Ay
    Mz=-M_Az
    Px=abs(-Rx/2)
    Py=abs(Mz/(h+t1)-Ry/2)
    Pz=abs(My/(h+t1)-Rz/2)
    return(Px,Py,Pz)

def Pullthrough(Fx,Fy,Fz,Mz,n,D2,e1,e2,e3,s2,t2,L):
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
    F_z = abs(F_z)
    F_y = abs(F_y)
    F_z = abs(F_z)
    F_y = abs(F_y)
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
    # print("K_bry: ", K_bry, " K_ty: ", K_ty, " K_t: ", K_t, " P_ty: ", P_ty, " P_y: ", P_y, " P_bry: ", P_bry, " F_y: ", F_y, " F_z: ", F_z)
    return SF




def BearingFailure(Ax, Az, My, D_2, t_2, L, n, sigmamaterial):
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

    return sigmamaterial / sigma


def Thermal(t_2, Ep, Eb, D2, Ap, Ab): # t1, E modulus lug, E modulus bolt, D2, Thermal coefficient lug, Thermal coefficient bolt,
    sa = 4* t_2 /(Ep*3.14*1.5*D2**2)
    sb = 1/ Eb * 4* t_2 /(3.14* D2**2/4)
    
    stupidletter = sa/(sa+sb)
    
    Force = (Ap - Ab)*250 * Eb *3.14 * D2**2/4 * (1-stupidletter)*(3.14*(D2/2)**2)
    
    return Force
    
    



def MassCalc(s2, D1, t1, w, t2, L, n, D2, rho):
    VolumeBP = t2*w*L - n*math.pi*((D2/2)**2)*t2
    VolumeF  = s2*w*t1 + 0.5*math.pi*((0.5*w)**2)*t1 - math.pi*((D1/2)**2)*t1

    Volume   = VolumeBP + 2 * VolumeF

    mass = Volume * rho

    print("VolumeBP, VolumeF, Volume, mass", VolumeBP, VolumeF, Volume, mass)

    return mass

# fortele=Forces(c.Fx, c.Fy, c.Fz, c.H, c.W)
# Fortax=fortele[12]
# Fortaz=fortele[13]
# My=fortele[10]

# LOOP FOR LUG A:
# Author: Seppe
# Enter general forces here: (Taken from the report)
Fx = 1.953*215.82
Fy = 1.953*215.82
Fz = 1.953*(-686.7)

# Then the forces at lug A are calculated:
Ax = Forces(Fx,Fy,Fz,0.450,0.975)[0]
Ay = Forces(Fx,Fy,Fz,0.450,0.975)[1]
Az = Forces(Fx,Fy,Fz,0.450,0.975)[2]

M_Ay = Forces(Fx,Fy,Fz,0.450,0.975)[6]
M_Az = Forces(Fx,Fy,Fz,0.450,0.975)[7]

# Then we also need some material properties:
Material = ["7075-T6", "2014-T6", "SEA-AISI 4340", "aged grade 250 maraging steel"]
Tau_max = [331*10**6, 290*10**6, 470*10**6, 1060*10**6] # https://asm.matweb.com/search/specificmaterial.asp?bassnum=ma7075t6, https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2014T6, https://www.azom.com/article.aspx?ArticleID=6772, https://www.makeitfrom.com/material-properties/Aged-Grade-250-Maraging-Steel
S_ty = [503*10**6, 414*10**6, 470*10**6, 1740*10**6]
rho = [2810, 2780, 7800, 8200]
E = [71.7*10*9, 73.1*10*9, 200*10*9, 190*10*9]
a = [23.6*10*(-6), 23.6*10*(-6), 12.3*10*(-6), 10.1*10*(-6)] #thermal expansion coefficient
ab = 8.6*10**(-6) # Thermal expansion of the bolt
Eb = 113.8 * 10 ** 9 # Young's modulus of the bolt 
# First iterate over flanges to determine s2

data = [['Iteration', 'D1','D2', 'L', 'W', 't1', 't2', 'n', 'SF Pullthrough', 'SF Flangefailure', 'SF Bearing', 'mass']]

iteration = 1

SFMAX = -1
massmax = 100
bestconfig = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


j = 3 # select material

for D2 in np.arange(0.002,0.003,0.0005):
    for L in np.arange((6*D2+4*0.001),(6*D2+4*0.001+0.05),0.01): # min 4x e1 + 4x t1
        for t2 in np.arange(0.0005,0.003,0.0005):
            for n in range(4,8,2):
                e1 = 1.5 * D2
                e3 = 2.5 * D2
                W = 2 * e1 + ((n/2)-1)* e3
                for D1 in np.arange(0.003,0.005,0.0005): 
                    TF = Thermal(t2, E[j], Eb, D2, a[j], ab) # extra "thermal force" on the fasteners
                    Pullthrougharray = Pullthrough(Ax, Ay+TF, Az, M_Az, n, D2, 1.5*D2, 1.5*D2, 2.5*D2, 2*D1, t2, L)
                    Tau_max_list = []
                    for i in range(int((n/2))):
                        Tau_max_list.append([Tau_max[j], Tau_max[j]])
                    Tau_max_array = np.array(Tau_max_list)
                    AbsPullthrougharray = abs(Pullthrougharray)
                    Marginsarray = Tau_max_array - AbsPullthrougharray
                    negative = np.any(Marginsarray<0)
                    if negative == False:
                        max = np.max(AbsPullthrougharray)
                        SFPullthrough = Tau_max[j] / max
                    else: 
                        SFPullthrough = 0.1
                
                    for t1 in np.arange(0.0005,0.0025,0.0005):
                        SFflange = FlangeFailure(W,D1,t1,S_ty[j],Ay,Az)
                        SFbearing = BearingFailure(Ax,Az,M_Ay,D2,t2,L,n,S_ty[j])

                        mass = MassCalc(2*D1, D1, t1, W, t2, L,n,D2, rho[j])

                        if mass <= massmax and SFbearing >=1.50 and SFflange >= 1.50 and SFPullthrough >= 1.50:
                            massmax = mass
                            bestconfig = [iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass]

                        data.append([iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass])
                        print("Iteration: ", iteration)
                        iteration += 1

# 7075:
# iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass
# 5917, 0.007, 0.004, 0.044, 0.022, 0.002, 0.001, 4, 2.1232695169359177, 1.837384353829659, 1.5211005492494596, 0.007744533938115955 - iterate t1
# 355, 0.004, 0.002, 0.02, 0.011, 0.005, 0.005, 4, 3.3568615327569638, 1.617837560704936, 1.7613556671596793, 0.006369343626288592 - t1 goes back, try 4
# 355, 0.004, 0.002, 0.028, 0.011, 0.005, 0.005, 4, 4.141821361462658, 1.617837560704936, 1.8374009395493285, 0.007605743626288591 - t1 stays high, try 4.5
# 43, 0.004, 0.002, 0.03, 0.011, 0.005000000000000001, 0.0045000000000000005, 4, 3.8843383663581728, 1.617837560704937, 1.6618407813805451, 0.007468849377001768 - yea that's fine

# 2014:
# iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass
# 356, 0.004, 0.002, 0.032, 0.011, 0.006, 0.005, 4, 3.8725886412641772, 1.5921552243557482, 1.5257219624476142, 0.008819746729987187 - iterate t1
# 356, 0.004, 0.002, 0.036000000000000004, 0.011, 0.006, 0.005, 4, 4.073277916556039, 1.5921552243557482, 1.5344719211000153, 0.009431346729987189 - works
# 327, 0.0045000000000000005, 0.0025, 0.039, 0.01375, 0.005, 0.004, 4, 4.369215930449027, 1.5549366235331576, 1.9365450632287333, 0.010806871244206805
# 683, 0.0045000000000000005, 0.0025, 0.035, 0.01375, 0.005000000000000001, 0.0035, 4, 3.662063344698583, 1.554936623533158, 1.6812942495977607, 0.009553626330384868 - final one

# SEA-AISI 4340
# iteration, D1, D2, L, W, t1, t2, n, SFPullthrough, SFflange, SFbearing, mass
# 362, 0.005, 0.002, 0.02, 0.011, 0.005, 0.005, 4, 4.590100897647908, 1.6114823300501229, 1.645799529950396, 0.018844679060487524 - t1 to 5
# 3019, 0.007, 0.003, 0.038000000000000006, 0.024, 0.002, 0.001, 6, 3.1414727734761145, 1.6755001765483954, 1.5771190909789377, 0.023122547318399085 - t1 jumps, try 4
# 362, 0.005, 0.002, 0.028, 0.011, 0.005, 0.005, 4, 5.603173917310993, 1.6114823300501229, 1.7168557486842633, 0.022276679060487522 - stays at config 1, will try 4.5
# 39, 0.0045000000000000005, 0.002, 0.03, 0.011, 0.005000000000000001, 0.0045000000000000005, 4, 5.37567477946074, 1.555156787724529, 1.5528134537750622, 0.021329677925422286 - stays at t1 = 5 anyway, will try 5
# 4, 0.0035, 0.002, 0.032, 0.011, 0.005500000000000001, 0.0045, 4, 5.844240547262905, 1.5162556680017683, 1.5588898311964752, 0.021772150978176922 - lol t1 increased, try 5.5
# 71, 0.0035, 0.002, 0.034, 0.011, 0.005499999999999999, 0.0045000000000000005, 4, 6.009301832350852, 1.5162556680017678, 1.5638015663042069, 0.02254435097817692 - use this one

# aged grade 250 maraging steel
# 57, 0.004, 0.002, 0.024, 0.016, 0.001, 0.001, 6, 3.508227425615884, 1.8439306335992622, 2.4707232695457675, 0.006536052987971815 - t1 to 1
# 57, 0.004, 0.002, 0.016, 0.016, 0.001, 0.001, 6, 2.2764795535768054, 1.8439306335992622, 2.1413070062304707, 0.005486452987971816
# 842, 0.004, 0.0025, 0.019, 0.01375, 0.001, 0.001, 4, 2.586284929809254, 1.809716427107034, 1.7764952488409, 0.0047967674886201395 - works



with open('Designpoints.csv', 'w', newline = '') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
    csvfile.close()

print(bestconfig)