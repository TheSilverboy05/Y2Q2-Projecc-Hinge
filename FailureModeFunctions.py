def Forces(Fx, Fy, Fz, H, W):
    """
    Returns a list [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz]
    Fx = Force in x-direction
    Fy = Force in y-direction
    Fz = Force in z-direction
    H  = Height of the solar panel in z-direction
    W  = length of the solar panel
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

    return [Ax, Ay, Az, Bx, By, Bz, M_Ay,M_Ay, M_Az, M_Bz]


def FlangeFailure(W,D,t,F_ty,F_y,F_z):

    #These need to be redone
    K_ty = 1 
    K_t = 1
    K_bru = 1
    
    A_br = D*t
    A_t = ((W-D)/t)

    P_ty = K_ty*A_br*F_ty
    R_tr = F_z/P_ty

    P_u = K_t*F_ty*A_t
    P_bru = K_bru*F_ty*A_br

    if P_u > P_bru:
        R_a = Fy/P_bru
    else:
        R_a = F_y/P_u

    SF = (1/((R_a^1.6+R_tr^1.6)^0.625))-1
    return SF


