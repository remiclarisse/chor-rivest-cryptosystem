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
    sieveSize = 4 * q ** 2
    nbIter = 0
    sieveTable = []
    constantsTable = []
    a, b, c, d = Fq2(0), Fq2(0), Fq2(0), Fq2(0)
    linear_poly = [ const for const in Fq2 ]
    while nbIter < sieveSize :
        while True :
            [a, b, c, d] = next_tuple ([a, b, c, d])
            print str((nbIter * 100 / sieveSize).n(digits=3)) + "\t" + str([a, b, c, d])
            if is_valid_quadruplet (a, b, c, d, q) :
                break
        aq = a ** q
        bq = b ** q
        cq = c ** q
        dq = d ** q
        P = Fq2X((c * aq - a * cq) * X * h0
               + (d * aq - b * cq) * h0
               + (c * bq - a * dq) * X * h1
               + (d * bq - b * dq) * h1)
        if is_split (P) :
            Q = Fq2X(h1 * (c * X + d))
            for gamma in range (q) :
                Q *= ((a - gamma * c) * X + b - gamma * d)
            indEqn = [ [linear_poly.index(poly.constant_coefficient()), mult] for poly, mult in list(factor(P)) ] + [ [linear_poly.index(poly.constant_coefficient()), -mult] for poly, mult in list(factor(Q)) ]
            indEqn.sort()
            nicen (indEqn)
            if indEqn not in sieveTable :
                sieveTable.append(indEqn)
                constantsTable.append(Q.lc() * P.lc() ** (-1))
                nbIter += 1
    return linear_poly, sieveTable, constantsTable

def nicen (L) :
    l = len (L)
    i = 0
    while i < l - 1 :
        if L[i][0] == L[i + 1][0] : # regrouper deux polynômes identiques
            L[i][1] += L[i + 1][1]
            L.pop(i + 1)
            l -= 1
            if L[i][1] == 0 : # si la multiplicité est nulle, le supprimer
                L.pop(i)
                i -= 1
                l -= 1
        i += 1

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
    order = card - 1
    fa = list(factor(order))
    for p, i in fa :
        order /= p ** i
        y = x ** order
        for j in range (i + 1) :
            if y == A.one() :
                break
            else :
                order *= p
                y = y ** p
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

def pohlig_hellman (h, g, n, fa) :
    a = [0 for i in range (len (fa))]
    for p, i in fa :
        g0 = g ** (n / p)
        ind = fa.index((p, i))
        for j in range (1, i + 1) :
            h0 = (  h * g ** (-a[ind])  ) ** (n / (p ** j))
            b = baby_step_giant_step (h0, g0, p)
            a[ind] = a[ind] + b * p ** (j - 1)
    moduli = [ p ** i for p, i in fa ]
    res = CRT (a, moduli)
    return res

def baby_step_giant_step (h, g, n) :
    m = int(ceil (sqrt (n)))
    L = []
    for i in range (m + 1) :
        u = g ** i
        L += [ hash(u) ]
    u = u ** (-1)
    y = h
    j = 0
    while hash(y) not in L :
        y = y * u
        j = j + 1
    return L.index(hash(y)) + m * j

def get_base_field_logs (g) :
    baseField = g.parent().base().base()
    order = baseField.cardinality() - 1
    fa = list(factor (order))
    orderg = multorder(g)
    g0 = g ** (orderg / order)
    logs = []
    for elem in baseField :
        if elem != baseField.zero() :
            l = pohlig_hellman (elem, g0, order, fa)
            logs += [ [ elem, mod (l * orderg / order, orderg) ] ]
    return zip(*logs)

def make_matrix (sieveTable, basisSize) :
    sieveSize = len(sieveTable)
    M = Matrix(ZZ, sieveSize, basisSize, sparse=True)
    for i in range (sieveSize) :
        for poly, mult in sieveTable[i] :
            M[i, poly] = mult
        print (i * 100 / sieveSize).n(digits=3)
    return M

def solve_logs_basis (M, constantsTable, baseFieldLogs, modulus) :
    B = [ baseFieldLogs[1][baseFieldLogs[0].index(elem)] for elem in constantsTable ]
    return M.change_ring(IntegerModRing(modulus)).solve_right(vector(B))

def joux_algorithm (p, h) :
    print "Embedding field..."
    Fp2h, h0, h1 = embed_field (p, h)
    print "Sieving linear polynomials..."
    linear_poly, sieveTable, constantsTable = sieving_linear_poly (p, h0, h1)
    print "Picking primitive element..."
    while True :
        g = Fp2h.random_element()
        if is_primitive(g) :
            break
    print "Computing logarithms of the base field..."
    baseFieldLogs = get_base_field_logs (g)
    print "Making equations out of the sieve..."
    make_matrix (sieveTable, len(linear_poly))
    print "Solving the system to recover the logarithms..."
    basisLogs = solve_logs_basis (M, constantsTable, baseFieldLogs, g.multiplicative_order())
    return baseFieldLogs, basis, basisLogs, g
