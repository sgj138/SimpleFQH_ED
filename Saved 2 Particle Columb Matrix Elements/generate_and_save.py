 
from math import comb, lgamma, gamma
import numpy as np

def vcolumb(l):
    return gamma(l+1/2)/gamma(l+1)



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


def create_2_particle_fermionic_basis(L):
    basis = []
    for i in range(L+1):
        if(i<L-i):
            state = [i,L-i]
            basis.append(state)
    return basis

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


V_blocks = {}
for l in range(99,150):
    basis = create_2_particle_fermionic_basis(l)
    dim = len(basis)
    V = np.zeros([dim,dim])
    for i,[m1_p,m2_p] in enumerate(basis):
        for j,[m1,m2] in enumerate(basis):
            V[i,j] = V_2(vcolumb,m1_p,m2_p,m1,m2)
    
    np.save(f"V2columb for L={l}.npy", V)

    # V_blocks[f"{l}"] = V
