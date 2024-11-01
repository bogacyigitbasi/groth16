## lets assuem we have a polynomial
## out = 4x^3 + 5x^2*y -7x*y^2 + 10x
import galois
import numpy as np
p = 71
GF = galois.GF(p)

# unit test values
x = GF(2)
y = GF(4)

out = 4*x**3 + 5*x**2*y - 7*x*y**2 + 10*x
print (out)

## we start with converting this polynomial to
## Arithmetic circuit to R1CS and then QAP

## R1CS
v1 = x*x
v2 = y*y
v3 = 5*v1*y
v4 = 4*x*v1
# out -v4 -v3 -10*x = -7*x*v2

# witness vector [1, out, x, y, v1, v2, v3, v4]
w = GF([1, out, x, y, v1, v2, v3, v4])

# L = [0,0,1,0,0,0,0,0]     R = [0,0,1,0,0,0,0,0]   O = [0,0,0,0,1,0,0,0]
#     [0,0,0,1,0,0,0,0]         [0,0,0,1,0,0,0,0]       [0,0,0,0,0,1,0,0]
#     [0,0,0,0,5,0,0,0]         [0,0,0,1,0,0,0,0]       [0,0,0,0,0,0,1,0]
#     [0,0,4,0,0,0,0,0]         [0,0,0,0,1,0,0,0]       [0,0,0,0,0,0,0,1]
#     [0,0,-7,0,0,0,0,0]        [0,0,0,0,0,1,0,0]       [0,1,-10,0,0,0,-1,-1]


# it must be L*R = O

L = GF([[0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,5,0,0,0],
        [0,0,4,0,0,0,0,0],
        [0,0,GF(p-7),0,0,0,0,0]])

print(L)

R = GF([[0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0]])
print(R)

O = GF([[0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,1],
        [0,1,GF(p-10),0,0,0,GF(p-1),GF(p-1)]])

print(O)

# Lw . Rw = Ow

Lw = np.dot(L, w)
Rw = np.dot(R, w)
Ow = np.dot(O, w)

mult = np.multiply(Lw, Rw)
print ("LwRw", mult)

print ("Ow", Ow)

assert np.all(mult == Ow)



