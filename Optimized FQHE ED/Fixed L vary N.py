
import numpy as np
from itertools import combinations
from math import gamma
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


def pairing_sign(p):
    inv = 0
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return 1 if inv % 2 == 0 else -1

def create_n_particle_fermionic_basis(n,L):
    basis = [state for state in list(combinations([i for i in range(L)],n)) if(sum(state)==L) ]
    
    return basis


def pairing(arr):
    n = len(arr)
    out = []
    for i in range(n):
        for j in range(i+1, n):
            pair = [arr[i], arr[j]]
            rest = [arr[k] for k in range(n) if k != i and k != j]
            out.append(pair + rest)
    return out
    
def vcolumb(l):
    return gamma(l+1/2)/gamma(l+1)

def vhardsphere(l):
    if(l<3):
        return 1
    else:
        return 0

def create_2_particle_fermionic_basis(L):
    basis = []
    for i in range(L+1):
        if(i<L-i):
            state = [i,L-i]
            basis.append(state)
    return basis


def indexof(m1,m2):
    basis = create_2_particle_fermionic_basis(m1+m2)
    for i in range(len(basis)):
        if(m1==basis[i][0] and m2==basis[i][1]):
            return i


N=3
m=3
L = int(0.5*m*N*(N-1))
V_blocks = {}

print(f"L={L}")

for l in range(L-2,L+3):
    basis = create_n_particle_fermionic_basis(N,l)
    dim = len(basis)
    V = np.zeros([dim,dim])
    for i,m in enumerate(basis):
        for j,m_p in enumerate(basis):
            sum_=0
            if(sum(m) == sum(m_p)):
                
                pairs = pairing(m)
                pairs_p = pairing(m_p)
        
                for pair in pairs:
                    for pair_p in pairs_p:
                        if(all(pair_p[i]==pair[i] for i in range(2,N))):
                            V_2 = np.load(f"Saved 2 Particle Matrix Elements\V2 for L={pair[0]+pair[1]}.npy")
                            sum_ += pairing_sign(pair)*pairing_sign(pair_p)*V_2[indexof(pair_p[0],pair_p[1]),indexof(pair[0],pair[1])]
                    
            
            V[i,j] = sum_
    V_blocks[f"{l}"] = V



    E = []
    for L_sector in V_blocks:
        # The martix blocks are sparse only for high block dimensions
        # if(int(L_sector) > 20):
        #     E_blocks,  eig_vec = eigsh(csr_matrix(V_blocks[L_sector]),k=10,which="SA")
        #     E.append(E_blocks)
        # else:
        #     E_blocks,  eig_vec = np.linalg.eigh(V_blocks[L_sector])
        #     E.append(E_blocks)
        
        E_blocks,  eig_vec = np.linalg.eigh(V_blocks[L_sector])
        E.append(E_blocks)


        # plt.plot([L_sector]*len(E_blocks), E_blocks+0.005*int(L_sector)**2, marker="_",markersize=10,ls="None", mew=2)
        plt.plot([L_sector]*len(E_blocks), E_blocks, marker="_",markersize=10,ls="None", mew=2)




plt.xlabel("Total Angular momentum L")
plt.ylabel("Eigenenergies E")

plt.title(str(N) +" Particle Columb Matrix Elements on a Disc | $L_{total}$=" + str(L))
# plt.legend()
plt.show()




  