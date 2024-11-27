import Configuration_I. as c
Fx=c.Fx
Fz=c.Fz
My=c.My


x=6   #this depends on the number of fasteners used and their positions
z=1



F_xbolt=Fx/n + My/(n*x)
F_zbolt=Fz/n + My/(n*z)



P = (F_xbolt**2+F_zbolt**2)**0.5

sigma=P/(D_2*t_2)    #max stress experienced by the bolt
#then compare sigma to the one of the material max strenght and see 

print(sigma/c.sigmamaterial)
