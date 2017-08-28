# p and h primes and h < p
def embed_field (p, h) :
    Fp2 = GF(p ** 2, 'a')
    FpY = PolynomialRing(GF(p), 'Y')
    Fp2X = PolynomialRing(Fp2, 'X')
    Ih, h0, h1 = pick_representation_polynome (p, h, FpY.gen())
    Fp2h = Fp2X.quotient_ring(Ih)
    return Fp2h, Fp2X(h0), Fp2X(h1)

def pick_representation_polynome (q, k, X) :
    A = X.parent()
    while true :
        h0 = A.random_element(2)
        h0 = h0 - h0(0)
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

def sieving_linear_poly (q, h0, h1) :
    X = h0.parent().gen()
    Fq2X = X.parent()
    Fq2 = Fq2X.base()
    sieveSize = q ** 3
    nbIter = 0
    hashTable = []
    sieveTable = []
    a, b, c, d = Fq2(0), Fq2(0), Fq2(0), Fq2(0)
    linear_poly = [ X + const for const in Fq2 ]
    while nbIter < sieveSize :
        while True :
            [a, b, c, d] = next_tuple ([a, b, c, d])
            if is_valid_quadruplet (a, b, c, d, q) :
                break
        print a, b, c, d
        P = Fq2X((c * a ** q - a * c ** q) * X * h0
               + (d * a ** q - b * c ** q) * h0
               + (c * b ** q - a * d ** q) * X * h1
               + (d * b ** q - b * d ** q) * h1)
        if is_split (P) :
            Q = Fq2X(h1 * (c * X + d))
            for gamma in range (q) :
                Q *= ((a - gamma * c) * X + b - gamma * d)
            indEqn = [ [linear_poly.index(poly), mult] for poly, mult in list(factor(P)) ] + [ [linear_poly.index(poly), -mult] for poly, mult in list(factor(Q)) ]
            indEqn.sort()
            print indEqn
            print list(factor(P)), list(factor(Q))
            if [hash(P), hash(Q)] not in hashTable :
                hashTable.append([hash(P), hash(Q)])
                sieveTable.append([ [(P.lc(), 1)] + list(factor(P)) ] + [ [(Q.lc(), 1)] + list(factor(Q)) ])
                nbIter += 1
                print (nbIter * 100 / sieveSize).n(digits=3)
    return sieveTable

def next_tuple (T) :
    K = T[0].parent()
    try :
        T[0] = K.next(T[0])
    except StopIteration :
        if len (T) == 1 :
            raise StopIteration
        return [K(0)] + next_tuple (T[1:])
    return T

def is_valid_quadruplet (a, b, c, d, p) :
    K = a.parent()
    if a * d == b * c :
        return False
    if a ** p == a and b ** p == b and c ** p == c and d ** p == d :
        return False
    return True

def multorder (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise ArithmeticError("ring not finite")
    return recursive_multorder (x, card - 1, prime_factors (card - 1))

def recursive_multorder (x, order, prime_fa) :
    A = x.parent()
    for p in prime_fa :
        if x ** (order / p) == A.one() :
            return recursive_multorder (x, order / p, [ f for f in prime_fa if (order / p) % f == 0 ])
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
    l = CRT (a, moduli)
    return h, l

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

def get_base_field_logs (g) :
    baseField = g.parent().base().base()
    multGpr = baseField.list()[1:]
    orderGpr = len(multGpr)
    orderExtMultGpr = multorder(g)
    g0 = g ** (orderExtMultGpr / orderGpr)
    fa = list(factor(orderGpr))
    res = list( pohlig_hellman ([ (g0, multGpr[i], fa) for i in range (orderGpr) ]) )
    res = [ ans for arg, ans in res ]
    logs = [ [ elem for elem, l in res ] ] + [ [ mod (l * orderExtMultGpr / orderGpr, orderExtMultGpr) for elem, l in res ] ]
    return logs

def make_equations (sieveTable, baseFieldLogs, modulus) :
    unknownsSide = []
    solutionSide = []
    for left, right in sieveTable :
        linear_poly = []
        sol = 0
        for poly, mult in left :
            if poly in baseFieldLogs[0] :
                sol -= baseFieldLogs[1][baseFieldLogs[0].index(poly)]
            else :
                linear_poly += [ (poly, mult) ]
        for poly, mult in right :
            if poly in baseFieldLogs[0] :
                sol += baseFieldLogs[1][baseFieldLogs[0].index(poly)]
            else :
                linear_poly += [ (poly, -mult) ]
        unknownsSide += [linear_poly]
        solutionSide += [sol]
    return unknownsSide, solutionSide

def make_basis (unknownsSide) :
    basis = []
    for line in unknownsSide :
        for poly, mult in line :
            basis.append(poly)
    basis = list(set(basis))
    return basis

def make_matrix (unknownsSide, basis, modulus) :
    M = Matrix(ZZ, len(unknownsSide), len(basis), sparse=True)
    sieveSize = len(unknownsSide)
    i = 0
    for line in unknownsSide :
        for poly, mult in line :
            M[i, basis.index(poly)] += mult
        i += 1
        print (i * 100 / sieveSize).n(digits=3)
    return M

def solve_logs_basis (M, solutionSide, modulus) :
    return M.change_ring(IntegerModRing(modulus)).solve_right(vector(solutionSide))

def joux_algorithm (p, h) :
    print "Embedding field..."
    Fp2h, h0, h1 = embed_field (p, h)
    print "Sieving linear polynomials..."
    sieveTable = sieving_linear_poly (p, h0, h1)
    print "Picking primitive element..."
    while True :
        g = Fp2h.random_element()
        if is_primitive(g) :
            break
    print "Computing logarithms of the base field..."
    baseFieldLogs = get_base_field_logs (g)
    print "Making equations out of the sieve..."
    modulus = p ** (2 * h) - 1
    unknownsSide, solutionSide = make_equations (sieveTable, baseFieldLogs, modulus)
    print "Making smoothness basis..."
    basis = make_basis (unknownsSide)
    print "Making the matrix..."
    M = make_matrix (unknownsSide, basis, modulus)
    print "Solving the system to recover the logarithms..."
    basisLogs = solve_logs_basis (M, solutionSide, modulus)
    return baseFieldLogs, basis, basisLogs, g
