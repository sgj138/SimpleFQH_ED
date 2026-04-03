
import numpy as np
from math import comb, lgamma, gamma
import matplotlib.pyplot as plt

def pairing_sign(p):
    inv = 0
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return 1 if inv % 2 == 0 else -1

def pairing(arr):
    n = len(arr)
    out = []
    for i in range(n):
        for j in range(i+1, n):
            pair = [arr[i], arr[j]]
            rest = [arr[k] for k in range(n) if k != i and k != j]
            out.append(pair + rest)
    return out


def create_2_particle_fermionic_basis(L):
    basis = []
    for i in range(L+1):
        if(i<L-i):
            state = [i,L-i]
            basis.append(state)
    return basis

def C(m1,m2,M,l):
    if(m1+m2==M+l):
        # print(np.sqrt((2**(m1+m2))))
        ll=lgamma(m1+1)+lgamma(m2+1)-lgamma(M+1)-lgamma(l+1)
        A = 1.0/ np.sqrt((2**(m1+m2) * np.exp(ll)))
        # print(m1,m2,A)
        sum_ = 0
        for k in range(max(0,l-m2),min(l,m1)+1):
            # for fermions we have two coeff substracted
            sum_+= comb(m1,k)*comb(m2,l-k)*(-1)**(l-k)
        return A*sum_
    else:
        return 0
    
def vcolumb(l):
    return gamma(l+1/2)/gamma(l+1)

def vhardsphere(l):
    if(l<3):
        return 1
    else:
        return 0

def V_2(v_l,m1_prime,m2_prime,m1,m2):
    """Returns the matrix elements for a pseudopotential v_l upto an angular momentum L"""
    
    # fermionic_basis_sorted_upto_L = create_fermionic_basis(L)
    sum_= 0
    if(m1+m2 == m1_prime+m2_prime):
        for l in range(m1+m2+1):
            M = m1+m2-l
            # since M+l=L therefore M = L-l
            fermionic_term = C(m1_prime,m2_prime,M,l)*C(m1,m2,M,l)- C(m1_prime,m2_prime,M,l)*C(m2,m1,M,l)- C(m2_prime,m1_prime,M,l)*C(m1,m2,M,l)+ C(m2_prime,m1_prime,M,l)*C(m2,m1,M,l)
            sum_+= fermionic_term*v_l(l)/(2)

    return sum_

  

def PlotColumb(L):

    V_blocks = {}
    for l in range(L+1):
        basis = create_2_particle_fermionic_basis(l)
        dim = len(basis)
        V = np.zeros([dim,dim])
        for i,[m1_p,m2_p] in enumerate(basis):
            for j,[m1,m2] in enumerate(basis):
                V[i,j] = V_2(vcolumb,m1_p,m2_p,m1,m2)
        V_blocks[f"{l}"] = V


    vl = [vcolumb(l) for l in range(1,L,2)]
    for j in range(len(vl)):
        # plt.plot([i for i in range(L)], [vl[j]]*L,color="black",alpha=0.1)
        plt.plot([L-1,L], [vl[j]]*2,color="black",alpha=0.4)
        plt.text(L,vl[j]-0.0025,"$\mathcal{v}_{"+ str(j+1)+"}$" )


    for L_sector in V_blocks:
        print(L_sector)
        E_blocks,  eig_vec = np.linalg.eigh(V_blocks[L_sector])

        plt.plot([L_sector]*len(E_blocks), E_blocks, marker="_",markersize=10,ls="None", mew=2)
        for i in range(len(E_blocks)+1):
            if(int(L_sector)%2==0):
                plt.text(int(L_sector)-0.25, E_blocks[i]-0.012,f"{2*i+2, int(L_sector) -2*i-1}",fontsize=6)
            if(int(L_sector)%2==1):
                plt.text(int(L_sector)-0.25, E_blocks[i]-0.012,f"{2*i+1, int(L_sector) -2*i}",fontsize=6)

        # plt.plot([L_sector]*len(E_blocks), E_blocks+0.1*int(L_sector)**2, marker="_",markersize=12,ls="None", mew=2)


    
    # plt.xlim(0,L+2)
    # plt.ylim(0,1.1)

    plt.xlabel("Total Angular momentum L")
    plt.ylabel("Eigenenergies E")

    plt.title("2 Particle Columb Matrix Elements on a Disc | $L_{total}$=" + str(L))
    # plt.legend()
    plt.show()




def LowLyingSpectra(L):
    # plt.xlim(0,L+2)
    # plt.ylim(0,1.1)

    plt.xlabel("Total Angular momentum L")
    plt.ylabel("Eigenenergies E")

    plt.title("2 Particle Columb Matrix Elements on a Disc | $L_{total}$=" + str(L))
    # plt.legend()
    plt.show()

PlotColumb(21)

# LowLyingSpectra(50)


