# p and n primes and n < p
def embedField (p, n) :
    Fp2 = GF(p ** 2, 'a')
    Fp2X = PolynomialRing(Fp2, 'X')
    Ik, h0, h1 = pickRepresentationPolynome (p, n, Fp2X.gen())
    Fp2n = Fp2X.quotient_ring(Ik)
    return Fp2n, h0, h1

def pickRepresentationPolynome (q, k, X) :
    A = X.parent()
    while true :
        h0 = A.random_element(2)
        h1 = A.random_element(2)
        fac = list((h1 * X ** q - h0).factor())
        deg = [ poly.degree() for poly, mult in fac ]
        if k in deg :
            break
    Ik = fac[deg.index(k)][0]
    return Ik, h0, h1
