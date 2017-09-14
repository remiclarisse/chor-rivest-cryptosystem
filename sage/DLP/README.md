# Some Algorithms to Solve DLP

Because Chor-Rivest cryptosystem uses discrete logarithms in finite fields, we have encountered several algorithms solving the discrete logarithm problem (a.k.a. DLP).

## Pohlig-Hellman's Algorithm

The file `pohlig-hellman.sage` contains an implementation of *Pohlig-Hellman's Algorithm*: Given a field GF(q) and **knowing the factorisation of q-1** (strong assumption), solve the DLP for each maximal sub-r-group of GF(q)\*, where r is a prime divisor of q-1, and reconstruct the original logarithm.

Here, the implementation uses the baby-step giant-step algorithm in the cyclic sub-r-groups.

### Example

In SageMath, after attaching or loading the file `pohlig-hellman.sage`, try:

    sage: K.<a> = GF(23 ** 7)
    sage: g = K.primitive_element()
    sage: h = K.random_element()
    sage: %time a = pohlig_hellman(h, g)
    sage: g ** a == h

The SageMath function computing discrete logarithms is `discrete_log(h, g)`. However, the function we implemented
can be used for an extension of an extension of a finite field (see below).

## Hellman-Reyneri's Algorithm

The file `hellman-reyneri.sage` contains an implementation of *Hellman-Reyneri's Algorithm*: Given a finite field GF(q^k), a generator `g` and a list of elements `h`, solve the DLP for each element in `h`.

In SageMath, the construction of the finite field GF(q^k) is tricky because `GF(q^k)` is viewed as GF(q)[X]/(P) where P is an irreducible polynomial of degree n (if q = p^m then n = m*k). We believe that this construction makes the computations remarkably slow, but it works.

### Example

In SageMath, after attaching or loading the file `hellman-reyneri.sage`, try:

    sage: q = 2 ** 3
    sage: k = 5
    sage: nb_tests = 10 ** 2
    sage: F.<a> = FiniteField (q)
    sage: R.<x> = PolynomialRing (F)
    sage: mu = R.irreducible_element (k)
    sage: K.<u> = R.quotient_ring (mu)
    sage: g = pick_primitive_element (K)
    sage: lo = [ randint (1, q ** k - 1) for i in range (nb_tests) ]
    sage: h = [ g ** i for i in lo ]
    sage: %time logs_hr = hellman_reyneri (g, h)
    sage: print "Did it work? " + str(lo == logs_hr)

### Remark

After the sieving phase, the system might not have a solution, depending on how many relevant equations where collected! That is something to be improved!!

## Joux's Algorithm

> Have a look at *Coppersmith's Algorithm* and the *Function Field Sieve* (FFS).

The file `joux-algorithm.sage` contains the begin of an implementation of *Antoine Joux's Algorithm*: it implements only what is necessary to determine the logarithm of linear polynomials.

### Example

To this date, the implementation we have made does not work...

---
## Bibliography

Stephen C. Pohlig and Martin E. Hellman, **An Improved Algorithm for Computing Logarithms over GF(p) and Its Cryptographic Signifiance**, in *IEEE Transactions on Information Theory*, 1978

Martin E. Hellman and Justin M. Reyneri, **Fast Computation of Discrete Logarithms in GF(q)**, in *Advances in Cryptology: Proceedings of CRYPTO '82*, 1982

Don Coppersmith, **Fast Evaluation of Logarithms in Fields of Characteristic Two**, in *IEEE transactions on information theory*, 1984

Leonard M. Adleman and Ming-Deh A. Huang, **Function Field Sieve Method for Discrete Logarithms over Finite Fields**, in *Information and Computation*, 1999

Antoine Joux and Reynald Lercier, **The Function Field Sieve in the Medium Prime Case**, in *Advances in Cryptology: EUROCRYPT 2006*, 2006

Antoine Joux, **A new index calculus algorithm with complexity L(1/4+o(1)) in small characteristic**, 2013
