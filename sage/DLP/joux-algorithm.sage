# p and h primes and h < p
def embed_field (p, n) :
    Fp2 = GF(p ** 2, 'a')
    Fp2X = PolynomialRing(Fp2, 'X')
    Ik, h0, h1 = pick_representation_polynome (p, n, Fp2X.gen())
    Fp2n = Fp2X.quotient_ring(Ik)
    return Fp2n, h0, h1

def pick_representation_polynome (q, k, X) :
    A = X.parent()
    while true :
        h0 = A.random_element(2)
        h1 = A.random_element(1)
        fa = list((h1 * X ** q - h0).factor())
        deg = [ poly.degree() for poly, mult in fa ]
        if k in deg :
            break
    Ik = fa[deg.index(k)][0]
    return Ik, h0, h1

def is_split (P) :
    fa = list(P.factor())
    deg = [ l[0].degree() for l in fa ]
    deg = set(deg)
    if deg == {1} :
        return True
    else :
        return False

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
            return recurcive_multorder (x, order / p, [ f for f in prime_fa if (order / p) % f == 0 ])
    return order

def is_primitive (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise ArithmeticError("ring not finite")
    card = card - 1
    for p in prime_factors(card) :
        if x ** (card / p) == A.one() :
            return False
    return True

@parallel
def pohlig_hellman (g, h, fa) :
    n = 1
    for p, i in fa :
        n = n * p ** i
    a = [0 for i in range (len (fa))]
    index = 0
    for p, i in fa :
        for j in range (1, i + 1) :
            g0 = g ** (n / (p ** j))
            h0 = (g0 ** (-a[index])) * (h ** (n / (p ** j)))
            if h0 != 1 :
                g0 = g ** (n / p)
                lg = baby_step_giant_step (g0, h0, p ** j)
                a[index] = a[index] + lg * p ** (j - 1)
        index += 1
    moduli = [ p ** i for p, i in fa ]
    res = CRT (a, moduli)
    return res

def baby_step_giant_step (g, h, n) :
    m = int(ceil (sqrt (n)))
    L = [ g ** i for i in range (m + 1) ]
    u = g ** (-m)
    y = h
    j = 0
    while y not in L :
        y = y * u
        j = j + 1
    return L.index(y) + m * j

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
