
import numpy as np
from itertools import product
import numpy as np
from math import comb, lgamma, gamma
import matplotlib.pyplot as plt
import sympy as sp

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

def create_fermionic_basis_3(L):
    basis = []
    for l in range(L+1):
        for i in range(0,l):
            for j in range(0,i):
                if(j<i-j<l-i):
                    state = [j,i-j,l-i]
                    basis.append(state)
                    print(state,sum(state))

    return basis

def create_fermionic_basis(n,L):
    """Returns a fermionic basis for n particles upto a total angular momentum L"""
    # Note that the constraint m_1<L and m_2<L has been used not m_1+m_2<L 
    # n= number of particles:
    # L= total angular momentum

    a = [i for i in range(0,L+1)]
    # take cartesian product to produce all possibilities of states for whole hilbert space
    res = list(product(a, repeat=n))
    # print(res)
    no_dup_basis = []
    for i in range(len(res)):
        # print(res[i])
        # print(f"after set: {set(res[i])}")7890-

        # set returns list after removing duplicates
        if(len(res[i])==len(set(res[i]))):
            no_dup_basis.append(res[i])
    

    L_values = [sum(l) for l in no_dup_basis]
    sorted_indices = np.argsort(L_values)
    no_dup_basis_sorted = [no_dup_basis[i] for i in sorted_indices]

    fermionic_basis = []
    for i in range(len(no_dup_basis_sorted)):
        is_true = 0
        for j in range(n-1):
            if(no_dup_basis_sorted[i][j]<no_dup_basis_sorted[i][j+1]):
                # print(i,no_dup_basis_sorted[i])
                is_true+=1
        if(is_true==n-1):
            fermionic_basis.append(no_dup_basis_sorted[i])
    fermionic_basis_up_to_L = []
    for i in fermionic_basis:
        if(sum(i)<=L):
            fermionic_basis_up_to_L.append(i)

    return fermionic_basis_up_to_L

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

def V_2(v_l,m1_prime,m2_prime,m1,m2):
    """Returns the matrix elements of two particles for a pseudopotential v_l upto an angular momentum L"""
    
    # fermionic_basis_sorted_upto_L = create_fermionic_basis(L)
    sum_= 0
    if(m1+m2 == m1_prime+m2_prime):
        for l in range(m1+m2+1):
            M = m1+m2-l
            # since M+l=L therefore M = L-l
            fermionic_term = C(m1_prime,m2_prime,M,l)*C(m1,m2,M,l)- C(m1_prime,m2_prime,M,l)*C(m2,m1,M,l)- C(m2_prime,m1_prime,M,l)*C(m1,m2,M,l)+ C(m2_prime,m1_prime,M,l)*C(m2,m1,M,l)
            sum_+= fermionic_term*v_l(l)/(2)

    return sum_


L=11
basis = create_fermionic_basis_3(L)
dim = len(basis)
V = np.zeros([dim,dim])

# print(basis)


# for i,[m1,m2,m3] in enumerate(basis):
#     print("m Basis: ", [m1,m2,m3])
#     for j,[m1_p,m2_p,m3_p] in enumerate(basis):
#         V_12 = 0
#         V_13 = 0
#         V_23 = 0
#         if(m1+m2+m3 == m1_p+m2_p+m3_p):
#             print("m prime Basis: ", [m1_p,m2_p,m3_p])
#             if(m1==m1_p):
#                 V_23 = pairing_sign([m1,m2,m3])*V_2(vcolumb,m2_p,m3_p,m2,m3)
#             if(m2==m2_p):
#                 V_13 = pairing_sign([m1,m3,m2])*V_2(vcolumb,m1_p,m3_p,m1,m3)
#             if(m3==m3_p):
#                 V_12 = pairing_sign([m2,m3,m1])*V_2(vcolumb,m1_p,m2_p,m1,m2)
            
#             print("V_12: ",V_12)
#             print("V_23: ",V_23)
#             print("V_13: ",V_13)

#             V[i,j] = V_12 + V_13 + V_23
#             print(V[i,j],"\n")

#     print("-------------------------------------------------------------")


for i,[m1,m2,m3] in enumerate(basis):
    print("m Basis: ", [m1,m2,m3])
    for j,[m1_p,m2_p,m3_p] in enumerate(basis):
        print("m_prime Basis: ", [m1_p,m2_p,m3_p])
        sum_ = 0
        if(m1+m2+m3 == m1_p+m2_p+m3_p):
            pairs = pairing([m1,m2,m3])
            pairs_p = pairing([m1_p,m2_p,m3_p])
            print("Pairs: ",pairs)
            print("Pairs Prime: ",pairs_p)
            print("Sign: ",[pairing_sign(i) for i in pairs])
            print("\n")
            for pair in pairs:
                for pair_p in pairs_p:
                    if(pair_p[2]==pair[2]):
                        print("INTERACTING Pair :",[pair[0],pair[1]])
                        print("INTERACTING Pair prime:",[pair_p[1],pair_p[0]])
                        print("\n")
                        sum_ += pairing_sign(pair)*pairing_sign(pair_p)*V_2(vcolumb,pair_p[0],pair_p[1],pair[0],pair[1])
            
            V[i,j] = sum_

    print("-------------------------------------------------------------")


# print(np.max(np.abs(V-V.T)))

# for i,[m1_p,m2_p] in enumerate(basis):
#     for j,[m1,m2] in enumerate(basis):
#         V[i,j] = V_2(vcolumb,m1_p,m2_p,m1,m2)


print(sp.latex(sp.Matrix(V)))




E, eig_vec = np.linalg.eigh(V) 
E_sorted = -np.sort(-E)



eapp=np.round(E_sorted,5)
# print(eapp)
energies=-np.sort(-np.unique(eapp))
e=list(energies)
eappl=list(eapp)

C=np.arange(len(energies))
for C,E in zip(C,energies):
    plt.text(C,E,str(eappl.count(E)))

plt.plot(energies, ".")
# plt.plot([vcolumb(2*i+1) for i in range(len(energies))],marker=".")


# plt.plot(E_sorted,".")
# print(sp.latex(sp.Matrix(V)))




plt.title(f"L={L}")
plt.grid()
plt.show()




