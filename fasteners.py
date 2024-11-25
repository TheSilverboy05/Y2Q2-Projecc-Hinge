import Configuration_I as c
Fx=c.Fx/2
Fy=c.Fy/2
F_xbolt=Fx/n
F_ybolt=Fy/n

P = (F_xbolt**2+F_ybolt**2)**0.5

sigma=P/(D_2*t_2)
#then compare sigma to the one of the material max strenght and see 

print(sigma)
