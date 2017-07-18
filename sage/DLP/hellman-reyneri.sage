def hellman_reyneri (q, k, F) :
    Ik = find_good_irreducible_poly (q, k, F)

def find_good_irreducible_poly (q, k, F) :
    R.<x> = PolynomialRing (F)
    P = [ F(0) for i in range (k) ]
    found = False
    while not found :
        P = next_poly (P)
        Q = R(x ** k + R(P))
        if Q.is_irreducible() :
            found = True
    return Q

def next_poly (P) :
    F = P[0].parent()
    if P[0] == 1 :
        [F(0)] + P[1:]
    else if F.next(P[0]) != 1 :
        return [F.next(P[0])] + P[1:]
    else :
        return [F(1)] + next_poly(P[1:])
