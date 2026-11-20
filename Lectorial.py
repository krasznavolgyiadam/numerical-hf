import numpy as np
from pylebedev import PyLebedev
import matplotlib.pyplot as plt

def main():
    N = (8,16,32,64,128)
    r,wr = gauss_chebychev(N)
    print(np.sum(wr*r**2*np.exp(-r))-2)

    leb = PyLebedev()
    u, wa = leb.get_points_and_weights(11)
    print(u.shape,wa.shape, wa.sum())

    print(4 * np.pi *np.sum(wa*u[:,0]**2))
    print(4 * np.pi /3)

    P = np.einsum("i,jk->ijk", r, u)
    x, y, z = P[..., 0,], P[..., 1], P[...,2]

    W = np.einsum('i,j->ij', wr * r**2, 4*np.pi * wa)

    print(integrate(W, np.exp(-(x**2+y**2+z**2))))
    print(np.pi**1.5)

    S_1s_1s = integrate(W, np.conj(psi_1s(x,y,z))*psi_1s(x,y,z))
    print(S_1s_1s)

    S_2s_2s = integrate(W, np.conj(psi_2s(x,y,z))*psi_2s(x,y,z))
    print(S_1s_1s)

    S_2px_2px = integrate(W, np.conj(psi_2px(x,y,z))*psi_2px(x,y,z))
    print(S_1s_1s)


    S_aa_1 = []
    for N in N:
        
        r,wr = gauss_chebychev(N)
        P = np.einsum("i,jk->ijk", r, u)
        W = np.einsum('i,j->ij', wr * r**2, 4*np.pi * wa)

        S_1s_1s = integrate(W, np.conj(psi_1s(x,y,z))*psi_1s(x,y,z))
        S_aa_1.append(np.absolute(S_1s_1s-1))

    plt.plot(N,)

def integrate(W,f):
    return np.einsum("ij,ij->",W, f)

def psi_1s(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    R = 2*np.exp(-r)
    Y = 1/(np.sqrt(4*np.pi))
    return R*Y

def psi_2s(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    R = 1/(2*np.sqrt(2))*(2-r)*np.exp(-r)
    Y = np.sqrt(3/(4*np.pi))*x/r
    return R*Y

def psi_2px(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    R = 1/(2*np.sqrt(6))*r*np.exp(-r/2)
    Y = np.sqrt(3/(4*np.pi))*x/r
    return R*Y


def gauss_chebychev(N, rm=1.0):
    """Points r and weights w such that int_0^inf g(r) dr ~ sum(w* g(r))."""
    z = np.arange(1, N+1)
    x = np.cos(np.pi * z / (N+1))
    r = rm * (1+x)/ (1-x)
    w = (np.pi / (N+1) * np.sin(np.pi*z/(N+1))**2*2*rm/(np.sqrt(1-x**2)*(1-x)**2))
    return r, w

if __name__ == "__main__": 
    main()