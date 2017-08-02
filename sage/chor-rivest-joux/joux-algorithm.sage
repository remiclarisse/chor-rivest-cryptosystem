# p and h primes and h < p
def joux_algorithm (g, h) :
    return 1

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
        h1 = A.random_element(1)
        fac = list((h1 * X ** q - h0).factor())
        deg = [ poly.degree() for poly, mult in fac ]
        if k in deg :
            break
    Ik = fac[deg.index(k)][0]
    return Ik, h0, h1

def sieving_linear_poly (q, X, h0, h1) :
    Fq2X = X.parent()
    Fq2 = Fq2X.base()
    sieveSize = q ** 2
    nbIter = 0
    hashTable = []
    sieveTable = []
    while nbIter < sieveSize :
        a = Fq2.random_element()
        b = Fq2.random_element()
        c = Fq2.random_element()
        d = Fq2.random_element()
        if a * d != b * c :
            P = Fq2X((c * a ** q - a * c ** q) * X * h0
                   + (d * a ** q - b * c ** q) * h0
                   + (c * b ** q - a * d ** q) * X * h1
                   + (d * b ** q - b * d ** q) * h1)
            if is_split (P) :
                Q = Fq2X(c * X + d)
                for gamma in range (q) :
                    Q *= ((a - gamma * c) * X + b - gamma * d)
                if [hash(P), hash(Q)] not in hashTable :
                    hashTable.append([hash(P), hash(Q)])
                    sieveTable.append([ [(P.lc(), 1)] + list(factor(P)) + [(h1.lc(), 1)] + list(factor(h1))] + [ [(Q.lc(), 1)] + list(factor(Q)) ])
                    nbIter += 1
                    print (nbIter * 100 / sieveSize).n(digits=3)
    return sieveTable

def get_unknowns (sieveTable) :
    unknowns = []
    i = 0
    sieveSize = len(sieveTable)
    for left, right in sieveTable :
        i += 1
        for poly, mult in left :
            unknowns.append(poly)
        for poly, mult in right :
            unknowns.append(poly)
        print (i * 100 / sieveSize).n(digits=3)
    unknowns = list(set(unknowns))
    return unknowns

def make_matrix_relation (sieveTable, basis) :
    sieveSize = len(sieveTable)
    M = Matrix(ZZ, sieveSize, len(basis), sparse=True)
    i = 0
    for P, Q in sieveTable :
        for poly, mult in P :
            M[i, basis.index(poly)] = mult
        for poly, mult in Q :
            M[i, basis.index(poly)] = -mult
        i += 1
        print (i * 100 / sieveSize).n(digits=3)
    return M

def solv_logs (M) :


def multorder (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise ArithmeticError("ring not finite")
    return recurcive_multorder (x, card - 1, prime_factors (card - 1))

def recurcive_multorder (x, order, prime_fa) :
    A = x.parent()
    for p in prime_fa :
        if x ** (order / p) == A.one() :
            return recurcive_multorder (x, order / p, [ fa for fa in prime_fa if (order / p) % fa == 0 ])
    return order

def is_split (P) :
    list_fact = list(P.factor())
    deg = [ l[0].degree() for l in list_fact ]
    deg = set(deg)
    if deg == {1} :
        return True
    else :
        return False
