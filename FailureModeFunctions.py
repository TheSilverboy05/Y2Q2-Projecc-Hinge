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

