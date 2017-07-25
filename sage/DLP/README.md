# Some Algorithms to Solve DLP

## Pohlig-Hellman

The file `pohlig-hellman.sage` contains an implementation of *Pohlig-Hellman Algorithm* : Given a field Z/pZ and **knowing the factorisation of p-1** (strong assumption), solve the DLP for each maximal sub-r-group of (Z/pZ)\*, where r is a divisor of p-1, and reconstruct the original logarithm.

Here, it uses an exhaustive seach in the cyclic sub-r-groups.

### Example

In `sage`, after attaching or loading the file `pohlig-hellman.sage`, try :

    sage: p = 2^15 * 3^6 * 5^4 + 1 # p is a prime number
    sage: g = GF(p)(4107883536) # g is primitive modulo p
    sage: h = GF(p)(572067100)
    sage: Fa = [[2,15],[3,6],[5,4]]
    sage: log (h,g) # is equal to 2074276358
    sage: %time pohlig_hellman (g, h, Fa) # must equal the value above

## Hellman-Reyneri

The file `hellman-reyneri.sage` contains an implementation of *Hellman-Reyneri Algorithm* : Given a finite field GF(q^k), a generator `g` and a list of elements `h`, solve the DLP for each element in `h`.

In `sage`, the construction of the finite field GF(q^k) is tricky because `GF(q^k)` is viewed as GF(p)[X]/(P) where P is an irreducible polynomial of degree n (if q = p^m then n = m*k). I believe that this construction makes the computation remarquably slow, but it works.

### Example

In `sage`, after attaching or loading the file `hellman-reyneri.sage`, try :

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

After the sieving phase, the system might not have a solution, depending on how many relevant equations where collected ! That is something to improve !!

---
## Bibliography

Stephen C. Pohlig and Martin E. Hellman, **An Improved Algorithm for Computing Logarithms over GF(p) and Its Cryptographic Signifiance**, in *IEEE Transactions on Information Theory*, 1978

Martin E. Hellman and Justin M. Reyneri, **Fast Computation of Discrete Logarithms in GF(q)**, in *Advances in Cryptology: Proceedings of CRYPTO '82*, 1982
