In this example we want to prove that we know 2 values x,y such that
4x^3 + 5x^2*y -7x*y^2 + 10x = (some value)

**Arithmetic Circuit**

- As a rule of thumb, protocol starts with the arithmetic circuit

**R1CS**

- Then we create R1CS form of the circuit, it has been proven that they are the same in the rareskills.io and I practiced in the repository in the past as well.

- When we have the R1CS matrices, witness vector we check the matrix multiplication such that
  Lw.Rw = Ow assert if not, cause we havent done anything yet. Its just matrix multiplication.

**Lagrange Interpolation**

- We simply divide each matrix columns' as vectors. This works thanks to homomorphism. Based on the number of constraints in the R1CS circuit we pick points on x axis to interpolate the points respectively.

- Please note that, Lagrange polynomial/interpolation is used to find the polynomial that intersects/interpolates the points given. So what we are doing with the previous step is to find the polynomial equivalents of the matrixes, in a way they are polynomials now and will be easier to compare the values we have.

**Schwartz Zippel Lemma**
One important step here is to introduce succinctness to the algorithm, Schwartz-Zippel lemma states that:

- Lets have two polynomials f with degree dx and g with degree dy, these two polynomials can intersect at max(dx, dy) points.
- If we evaluate a random point(tau) on these two polynomials, if the f(tau) == g(tau) is highly likely f=g (considering the degrees dx, dy << Finite P)

What it means for groth16 is, we dont have to compare all one by one, nxm multiplied matrices. We have to calculate a polynomial value once. and verify its correctnes in O(1) time not O(n) (or nxm to be more accurate.)

**Elliptic Curve**

- An EC is a 3rd order polynomial, we use them to generate keys and increase the security of the systems. On an EC, if you pick a number and multiply it with a scalar you get another point on the curve. It doesnt matter how many times you have done that, you end up with a new point. However, it is impossible due to discrete logarithm to find the value n.
- Adding two points yield another one on the curve.
- There are some special 3rd order polynomials that support bilinear pairing, we have to use them on the protocol to enable some amazing features.

**Trusted Setup**

- Trusted setup provides required global parameters, simply its picking a random value tau and computing a reference string to it by multiplying with a point on a curve for both prover and the verifier.
- Once its done we have a list of points multiplied by scalar (tau) both parties can have it and operate.
- Prover will use it to generate the proof without having the value of tau, verifier will use it to make sure that prover is indeed proving the correct thing. This is a risk for verifier, if he has no evidence of what is being proven, then a prover can prove an arbitrary equation to scam.
