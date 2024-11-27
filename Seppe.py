import numpy as np

def BoltsLoad(Fx,Fy,Fz,Mz,n,D2,e1,e2,e3,s2,t2,W,L):
    # Author: Seppe
    rows = int(n/2)
    columns = 2
    
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
        # output in array
    # Force closest to the neutral axis:
    denominatorfactor = 0
    Fzlist = []

    if (n%2) == 0:
        for i in range(int(n/4)):
            denominatorfactor += 2*i+1
        Fclosest = ((Fz*(s2+(t2/2)))/(2*e3*denominatorfactor))

        for i in range(int(n/4+1),int(-(n/4+2)),-2):
            Fzlist.append([Fclosest*i, Fclosest*i])
    

    elif (n%2) != 0:
        print("ik stink")

    Fzarray = np.array(Fzlist)

    # Forces due to Mz
        # output in array

    # Sum of all forces per bolt
    Farray = Fxarray + Fyarray
        # Add all arrays
    
    return(Fzarray)


print(BoltsLoad(100,100,100,100,8,10,30,30,30,50,5,100,200))
