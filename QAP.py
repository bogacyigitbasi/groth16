## lets assuem we have a polynomial
## out = 4x^3 + 5x^2*y -7x*y^2 + 10x
# and we know 2 secret values that satisfies the equation.
# x=2, y=4, result=49
import numpy as np

import galois

p = 79
GF = galois.GF(p)
# unit test values
x = 2
y = 4

# out = 4*x**3 + 5*x**2*y - 7*x*y**2 + 10*x
# print (out)

## we start with converting this polynomial to
## Arithmetic circuit to R1CS and then QAP

## R1CS
v1 = x*x
v2 = y*y
v3 = 5*v1*y
v4 = 4*x*v1
out = -7*x*v2 +v4 +v3 +10*x

# witness vector [1, out, x, y, v1, v2, v3, v4]
w = np.array([1, out, x, y, v1, v2, v3, v4])

# L = [0,0,1,0,0,0,0,0]     R = [0,0,1,0,0,0,0,0]   O = [0,0,0,0,1,0,0,0]
#     [0,0,0,1,0,0,0,0]         [0,0,0,1,0,0,0,0]       [0,0,0,0,0,1,0,0]
#     [0,0,0,0,5,0,0,0]         [0,0,0,1,0,0,0,0]       [0,0,0,0,0,0,1,0]
#     [0,0,4,0,0,0,0,0]         [0,0,0,0,1,0,0,0]       [0,0,0,0,0,0,0,1]
#     [0,0,-7,0,0,1,0,0]        [0,0,0,0,0,1,0,0]       [0,1,-10,0,0,0,-1,-1]


# it must be L*R = O

L = np.array([[0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,5,0,0,0],
        [0,0,4,0,0,0,0,0],
        [0,0,-7,0,0,0,0,0]])

# print(L)

R = np.array([[0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0]])
# print(R)

O = np.array([[0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,1],
        [0,1,-10,0,0,0,-1,-1]])

# print(O)

# Lw . Rw = Ow
assert all (np.equal(np.matmul(L, w)*np.matmul(R, w), np.matmul(O, w)))
# print ("LwRw", mult)

# print ("Ow", Ow)

# assert np.all(mult == Ow)

# QAP -> motivation here is to convert from vector field to polynomial field
# then we will apply schwartz zippel lemma to evaluate polynomials are identical
# L.R = 0 is becoming l(x).r(x) = o(x)
# Lagrange interpolation is used to find the polynomial that interpolates given points.
# we have 5 constraints/rows therefore we will interpolate on 5 points x = 1,2,3,4,5
# and we use each column as a vector of coefficients for each matrix.
print("-"*10)
L_GF = GF((L+p)%p)
# print(L_GF)

R_GF = GF((R+p)%p)
# print(R_GF)

O_GF = GF((O+p)%p)
# print(O_GF)
print("-"*10)

v1 = x*x
v2 = y*y
v3 = GF(5)*v1*y
v4 = GF(4)*x*v1
out = GF(-7+p)*x*v2 +v4 +v3 +GF(10)*x

# witness vector [1, out, x, y, v1, v2, v3, v4]
w = GF(np.array([1, out, x, y, v1, v2, v3, v4]))

def interpolate_columns(col):
    xs = GF(np.array([1,2,3,4,5]))
    mat = galois.lagrange_poly(xs, col)
    # print (mat)
    return mat
print("-"*10)
U_polys = np.apply_along_axis(interpolate_columns, 0, L_GF)
print("U: " ,U_polys)
V_polys = np.apply_along_axis(interpolate_columns, 0, R_GF)
print("V: " ,V_polys)
W_polys = np.apply_along_axis(interpolate_columns, 0, O_GF)
print("W: " ,W_polys)
print("-"*10)

assert all(np.equal(np.matmul(L_GF, w)*np.matmul(R_GF, w), np.matmul(O_GF, w)))


# tau = GF(20)
# print("U(tau):", [poly(tau) for poly in U_polys])
# print("V(tau):", [poly(tau) for poly in V_polys])
# print("W(tau):", [poly(tau) for poly in W_polys])

# assert U_polys(tau)*V_polys(tau) == W_polys(tau)

def inner_product_polys_witness(polys, witness):
    total = GF(0)
    for p,w in zip (polys, witness):
        total += p*w
    return total

print("-"*10)
Ux= inner_product_polys_witness(U_polys, w)
print("Ux: ",Ux)
Vx= inner_product_polys_witness(V_polys, w)
print("Vx: ", Vx)
Wx= inner_product_polys_witness(W_polys, w)
print("Wx: ",Wx)
# assert all(np.equal(np.matmul(L_GF, w)*np.matmul(R_GF, w), np.matmul(O_GF, w)))
print("-"*10)



# trusted setup needs to be managed by a trusted agent who generates tau
# and deletes after computing G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
# we will use bn128

from py_ecc.bls12_381 import G1, G2, multiply, add, pairing, eq, curve_order, is_on_curve, b,b2

from functools import reduce

### QAP Balancing T(x)
# it fails cause an imbalance arises on the multiplication of L and R and the
# equation becomes skewed. Left has a polynomial with degree 8, right with 4
# add the polynomial T(x)*H(x) where T(x)= (x -1) (x-2)(x-3)(x-4)(x-5)interpolated points

Tx = galois.Poly([1, p-1], field=GF)
for i in range(2,6):
    Tx *= galois.Poly([1,p-i],field=GF)
    # print("Tx1: ", Tx)

print("Tx: ",Tx)
# print(Tx(11))
tau = GF(20) # must be less than 79
# # T_tau = Tx(tau)
# # verify trusted setup
# tau_int = int(tau) % curve_order # the pairing function and elliptic curve operations in the py_ecc library require standard integers, not elements from galois.GF
# #
# pairing_left = pairing(multiply(G2, tau_int), multiply(G1, tau_int**2))
# pairing_right = pairing(multiply(G2,tau_int), multiply(G1,tau_int))**2

# assert pairing_left == pairing_right, "verification failed, ceremony is not trusted"

tau_int = int(tau)
# powers of tau for A , C on G1 # G1[τ^0], G1[τ^1], ..., G1[τ^d-1]
srs_tau_G1 = [multiply(G1, (tau_int**i)) for i in range (5)]
# print("PoTG1", srs_tau_G1)

# powers of tau for B on G2
# G2[τ^0], G2[τ^1], ..., G2[τ^d]
srs_tau_G2 = [multiply(G2, (tau_int**i)) for i in range (5)]
# print("PoTauG2: ",srs_tau_G2)

# G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
srs_tau_hx_tx_G1 = [multiply(G1, int(tau_int**i*(Tx(tau)))) for i in range (4)]
# print("PoT HxTx: ",srs_tau_hx_tx_G1)

print("srs_tau_G1 points:", srs_tau_G1)
print("srs_tau_G2 points:", srs_tau_G2)
print("srs_tau_hx_tx_G1 points:", srs_tau_hx_tx_G1)


# now lets check if
# assert Ux(tau)*Vx(tau) == Wx(tau) # fails as expected
## QAP Balancing H(x) and T(x)
# now we need the remaining polynomial H(x)
# U(x) * V(x) = W(x) + T(x)+H(x)
# H(x) = U(x) * V(x) - W(x) / T(x)
print("-"*10)

Hx = (Ux * Vx - Wx) // Tx

assert Ux(tau)*Vx(tau) == Wx(tau) + Tx(tau)*Hx(tau)
remainder = (Ux * Vx - Wx) % Tx

print("Hx degree: ",Hx.degree)
print("Tx degree: ",Tx.degree)
print("remainder: ", remainder)
print("-"*10)
# compute tau
print("Hx: ",Hx(tau))
print("-"*10)
### Prover send 3 EC points [A]G1, [B]G2, [C]G1
print("Evaluations at tau for Ux, Vx, Wx, Tx, and Hx:", Ux(tau), Vx(tau), Wx(tau), Tx(tau), Hx(tau))

print("-"*10)


def polynomial_hiding(poly, PoTau, G):
    coefficients = poly.coefficients()
    print("coefficients: ", coefficients)
    assert len(coefficients) == len(PoTau), "Length mismatch"
    result = multiply(G, 0)
    for i in range(len(PoTau)):
        result = add(result, multiply(PoTau[i], int(coefficients[-(i+1)])))
    return result

# [A]G1
A_G1 = polynomial_hiding(Ux, srs_tau_G1, G1)

# [B]G2
B_G2 = polynomial_hiding(Vx, srs_tau_G2, G2)

# [C]G1
W_G1 = polynomial_hiding(Wx, srs_tau_G1, G1)


HT_G1 = polynomial_hiding(Hx, srs_tau_hx_tx_G1, G1)
C_G1 = add(W_G1, HT_G1)
proof = [A_G1, B_G2, C_G1]

# print("Proof: ",proof)
print("-"*10)

assert is_on_curve(proof[0], b)
assert is_on_curve(proof[1], b2)
assert is_on_curve(proof[2], b)

#### Verifier will check the pairing
# it must be the same e(A, B) == e (C, G2)
# result = pairing(proof[1], proof[0]) == pairing(G2,proof[2])
# print("e(A, B) == e(C, G2[1]):", result)
# print("The proof is valid?", result)

pairing_left = pairing(proof[1], proof[0])
pairing_right = pairing(G2, proof[2])
print("Left pairing (e(A, B)):", pairing_left)
print("Right pairing (e(C, G2)):", pairing_right)
print("Pairing Equality Check:", pairing_left == pairing_right)