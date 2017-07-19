# attach("~/chor-rivest-cryptosystem/sage/DLP/hellman-reyneri.sage")
# q = 2 ** 2
# k = 5
# F.<a> = FiniteField (q)
# R.<x> = PolynomialRing (F)
# mu = R(x^5 + a*x^3 + a); mu.is_irreducible()
# K.<u> = R.quotient_ring (mu)
# hellman_reyneri (q, k, F, None, None, mu)

def hellman_reyneri (q, k, F, g, h, mu) :
    b = min (k, 10)
    R = PolynomialRing (F, 'x')
    Ik = find_good_irreducible_poly_S (q, k, b, F, R)
    K = R.quotient_ring(mu, 'u')
    Q = R.quotient_ring(Ik, 'v')
    u = K.gen()
    v = Q.gen()
    i = 0
    while (v ** i).minpoly() != mu and i < q**k :
        i += 1
    if i == q ** k :
        return "fail"
    t = v ** i
    M = Matrix( [ list(t**i) for i in range (k) ] ).transpose() # matrice de t à u (exprime t dans la nase 1, u, u^2, ..., u^{k-1})
    return K(t), M, Ik

def find_good_irreducible_poly_S (q, k, b, F, R) :
    x = R.gen()
    P = [ F(0), F(0), F.gen() ] + [ F(0) for i in range (b - 3) ]
    found = False
    while not found :
        P = next_poly (P)
        Q = R(x ** k + R(P))
        if Q.is_irreducible() :
            return Q
    return None

def next_poly (P) :
    F = P[0].parent()
    if P[0] == 1 :
        return [F(0)] + P[1:]
    elif F.next(P[0]) != 1 :
        return [F.next(P[0])] + P[1:]
    else :
        return [F(1)] + next_poly(P[1:])
