import numpy as np

def BM4_TH(x,K0,Kp,Kdp,V0,P0,T0,a,b,c):
    vol, Temp = x
    a1 = 3*V0*P0
    a2 = 3*V0*(3*K0-5*P0)/2
    a3 = V0*(9*K0*Kp - 36*K0 + 35*P0)/2
    a4 = 3*V0*(9*(K0**2)*Kdp + 9*K0*(Kp**2) - 63*K0*Kp + 143*K0 - 105*P0)/8.0
    f = 1/2*((V0/vol)**(2/3)-1)
    dfdV = - (V0**(2/3))/3/(vol**(5/3))
    Pc = (a1 + 2*a2*f + 3*a3*(f**2) + 4*a4*(f**3))*(-dfdV)
    Pth = (a - b*(vol/V0)+c*((vol/V0 )**2))/1000. *(Temp-T0)
    return Pc+Pth
