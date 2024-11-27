import math


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
    Az = -Fz/2
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

    return [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz, My, Mz]


def FlangeFailure(W,D,t,S_ty,F_y,F_z):
    A_br = D*t
    A_t = (W-D)*t

    # K_bry, fig. D1.14 in Bruhn, taken max curve with out of bounds if the relevant curve diverges
    if (0.5*W)/D < 3.3*(D/t)+0.57:
       K_bry = -1.25*(D/t)-0.998*(D/t)^2+0.11*(D/t)^3
    else:
        print("e/D conbined with D/t is out of bounds")
        K_bry = 100000000

    # K_ty, fig. D1.15 in Bruhn, only curve 3
    A_2 = 0.5*(W-D)
    A_1 = A_2+0.5*D-0.5*D*math.cos(45*(math.pi/180))
    A_av = 6/((a/A_1)+(2/A_2))
    K_ty = -4.72*10^(-3)+1/39(A_av/A_br)-0.341(A_av/A_br)^2

    # K_t, fig. D1.12 in Bruhn, only linear part otherwise out of bounds
    if W/D <= 2.9: 
        K_t = (W/D)*(-0.08/1.4)+1.0857
    else:
        K_t = 100000000
        print("W/D out of bounds (over 2.9)")

    # transverse
    P_ty = K_ty*A_br*S_ty
    R_tr = F_z/P_ty

    # axial
    P_y = K_t*S_ty*A_t
    P_bry = K_bry*S_ty*A_br

    if P_y > P_bry:
        R_a = F_y/P_bry
    else:
        R_a = F_y/P_y

    SF = (1/((R_a^1.6+R_tr^1.6)^0.625))-1
    return SF




def BearingFailure( D_2, t_2, n)
    
    F_xbolt=Fx/n + My/(n*x)
    F_zbolt=Fz/n + My/(n*z)

    P = (F_xbolt**2+F_zbolt**2)**0.5

    sigma=P/(D_2*t_2)    #max stress experienced by the bolt
    #then compare sigma to the one of the material max strenght and see how to lighten up the hinge

    sigmamaterial=1 #to be completed

    return sigma/sigmamaterial

