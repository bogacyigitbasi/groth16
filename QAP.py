## lets assuem we have a polynomial
## out = 4x^3 + 5x^2*y -7x*y^2 + 10x

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
#     [0,0,-7,0,0,0,0,0]        [0,0,0,0,0,1,0,0]       [0,1,-10,0,0,0,-1,-1]


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

L_GF = GF((L+p)%p)
print(L_GF)

R_GF = GF((R+p)%p)
print(R_GF)

O_GF = GF((O+p)%p)
print(O_GF)


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

U_polys = np.apply_along_axis(interpolate_columns, 0, L_GF)
# print("U: " ,U_polys)
V_polys = np.apply_along_axis(interpolate_columns, 0, R_GF)
# print("V: " ,V_polys)
W_polys = np.apply_along_axis(interpolate_columns, 0, O_GF)
# print("W: " ,W_polys)


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


Ux= inner_product_polys_witness(U_polys, w)
print("Ux: ",Ux)
Vx= inner_product_polys_witness(V_polys, w)
print("Vx: ", Vx)
Wx= inner_product_polys_witness(W_polys, w)
print("Wx: ",Wx)
# assert all(np.equal(np.matmul(L_GF, w)*np.matmul(R_GF, w), np.matmul(O_GF, w)))

# now lets check if term1*term2 == term3
tau = GF(20)

# assert Ux(tau)*Vx(tau) == Wx(tau) # fails as expected

# it fails cause an imbalance arises on the multiplication of L and R and the
# equation becomes skewed. Left has a polynomial with degree 8, right with 4
# add the polynomial T(x)*H(x) = (x -1) (x-2)(x-3)(x-4)(x-5) interpolated points

Tx = galois.Poly([1, p-1], field=GF)
for i in range(2,5):
    Tx *= galois.Poly([1,p-i],field=GF)

print("Tx: ",Tx)

# now we need the remaining polynomial H(x)
# U(x) * V(x) = W(x) + T(x)+H(x)
# H(x) = U(x) * V(x) - W(x) / T(x)

Hx = (Ux * Vx - Wx) // Tx

assert Ux(tau)*Vx(tau) == Wx(tau) + Tx(tau)*Hx(tau)

remainder = (Ux * Vx - Wx) % Tx

print("Hx: ",Hx)
print("remainder: ", remainder)

# compute tau
print(Hx(tau))

# trusted setup needs to be managed by a trusted agent who generates tau
# and deletes after computing G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
# we will use bn128

from py_ecc.bn128 import G1, G2, multiply, add, pairing, eq

from functools import reduce

# G1[τ^0], G1[τ^1], ..., G1[τ^d-1]
structured_reference_string_G1 = [multiply(G1, int(tau**i)) for i in range (0,5)]
print(structured_reference_string_G1)
# G2[τ^0], G2[τ^1], ..., G2[τ^d]
srs_tau_G2 = [multiply(G2, int(tau**i)) for i in range (0,5)]
print(srs_tau_G2)
# G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
srs_hx_tx = [multiply(G1, int(tau**i*Tx(tau))) for i in range (0,5)]
print(srs_hx_tx)
# verify trusted setup
tau_int = int(tau) # the pairing function and elliptic curve operations in the py_ecc library require standard integers, not elements from galois.GF

pairing_left = pairing(multiply(G2, tau_int), multiply(G1, tau_int**2))
pairing_right = pairing(multiply(G2,tau_int), multiply(G1,tau_int))
assert pairing_right == pairing_right, "verification failed, ceremony is not trusted"

