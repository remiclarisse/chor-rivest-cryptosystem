# reset()
# attach("/home/grace/rclariss/chor-rivest-cryptosystem/sage/chor-rivest.sage")
# attach("/home/grace/rclariss/chor-rivest-cryptosystem/sage/draft.sage")
# p = 197
# h = 24
# r = 8
# [PubKey, PrivKey] = CRGenerateKeys (p, h)
# [c, p, h, Q, alpha] = PubKey
# [t, g, sInv, d] = PrivKey
# gpr = g**((p**h-1)/(p**r-1))
# s = [sInv.index(i) for i in range (p)]


def blop (t, p, h) :
    K = t.parent()
    for j in range (h) :
        for a in range (1, p) :
            for b in range (p) :
                if t**(p**j) + a*t + b == 0:
                    print j, a, b
    return None

def blip (c, p, h, alpha, t) :
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [ int(mod (log (t + alpha[i], gg), r)) for i in range (p) ]
    C = [ int(mod (c[i]-c[0], r)) for i in range(p) ]
    A = [ int(mod (a[0]-a[i], r)) for i in range(p) ]
    for i in range (p) :
        d = gcd (C[i], r)
        for j in range (p) :
            if A[j] % d == 0 and C[i] != 0 and A[j] != 0:
                gamma = Integer(C[i] / d)
                modulus = Integer(r / d)
                beta = Integer(A[j] / d)
                L = ((gamma).inverse_mod(modulus) * (beta % modulus)) % modulus
                # print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d)
                for w in range (d) :
                    # if w % 100 == 0 and w != 0 :
                    #     print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d) + "\tw:" + str(w)
                    # if ok : sigma(0) = j and sigma(i) = 0
                    D = [ (L * C[l]) % r for l in range(p) ]
                    E = [ (a[k] - a[j]) % r for k in range(p) ]
                    if set(D) == set(E) :
                        return gg**L, int(mod(c[0] - log(t + alpha[j], gg**L), r - 1))
                    L = L + modulus
    return "fail"

def blup (c, p, h, alpha, pi, gpr, r) :
    n = h / r
    K = gpr.parent()
    A.<X> = PolynomialRing(K)
    Q = A(0)
    for i in range (n + 1) :
        L = A(1)
        for k in range (n + 1) :
            if k != i:
                L = L * A((alpha[pi[i]] - alpha[pi[k]]) ** (-1) * (X - alpha[pi[k]]))
        Q = Q + gpr ** c[i] * L
    R = Q.roots()
    T = []
    for root in R :
        T.append(-root[0])
    return T

def next_value (array, i, p) :
    array[i] = array[i] + 1
    while array[i] in array[:i] :
        array[i] = array[i] + 1
    if array[i] == p :
        next_value (array, i - 1, p)
        restart_value (array, i, p)

def restart_value (array, i, p) :
    array[i] = -1
    new_value = 0
    while new_value in array :
        new_value = new_value + 1
    array[i] = new_value

def increment (array, stone, p) :
    array[-1] = array[-1] + 1
    while array[-1] in array[:-1] :
        array[-1] = array[-1] + 1
    if array[-1] == p :
        next_value (array, -2, p)
        restart_value (array, -1, p)
    if array[1] != stone :
        return False
    else :
        return True

def blap (c, p, h, alpha, gpr, r, data) : # data === array of size two (e.g. [0, 1])
    n = h / r
    K = gpr.parent()
    sig = [-1 for i in range (p)] # pas obligé de le définir ici ?
    defined = [False for i in range (p)] # pas obligé de le définir ici ?
    newnbs = data + [i for i in range (2, n)] + [n - 1]
    found = False
    ok = True
    while not found and ok:
        ok = increment (newnbs, data[1], p)
        # print newnbs # verbose mode
        sig = newnbs + [-1 for i in range (p - (n + 1))]
        defined = [False for i in range (p)]
        for j in sig :
            if j != -1 :
                defined[j] = True
        carry_on = True
        l = 0
        while l < p and carry_on :
            if not defined[l] :
                Q = K(0)
                for j in range (n + 1) :
                    L = K(1)
                    for k in range (n + 1) :
                        if k != j :
                            L = L * (alpha[sig[j]] - alpha[sig[k]]) ** (-1) * (alpha[l] - alpha[sig[k]])
                    Q = Q + gpr ** c[j] * L
                i = 0
                while i < p and gpr ** c[i] != Q :
                    i = i + 1
                if i == p : # ordre des tests avec 'or' ?
                    carry_on = False
                elif sig[i] != -1 :
                    carry_on = False
                else :
                    sig[i] = l
                    defined[l] = True
            l = l + 1
        if carry_on :
            found = True
    if not ok :
        return "fail"
    else :
        return sig

def blipblip (c, p, h, alpha, t, pi) :
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [ int(mod (log (t + alpha[pi[i]], gg), r)) for i in range (p) ]
    C = [ Integer(mod (c[i]-c[0], r)) for i in range(p) ]
    A = [ int(mod (a[i]-a[0], r)) for i in range(p) ]
    for i in range (p) :
        if gcd (C[i], r) == 1 :
            L = (C[i]).inverse_mod(r) * (A[i])
            D = [ (L * C[l]) % r for l in range(p) ]
            if set(D) == set(A) :
                return gg**L, int(mod(c[0] - log(t + alpha[pi[0]], gg**L), r - 1))
    return "fail"

def crackMessage (PubKey, e, gpr, r) :
    # Récuperer la clé publique
    c = PubKey[0]
    p = PubKey[1]
    h = PubKey[2]
    q = p ** h
    alpha = PubKey[4]
    # Générer une clé privée équivalente
    pi = blap (c, p, h, alpha, gpr, r, [0, 1])
    piInv = [pi.index(i) for i in range (p)]
    T = blup (c, p, h, alpha, pi, gpr, r)
    t = T[0]
    mu = t.minimal_polynomial()
    g, d = blipblip (c, p, h, alpha, t, pi)
    # Faire le changement de base
    V = t.parent().vector_space()
    Minv = Matrix (GF(p), [V(t ** i) for i in range (h)]).transpose().inverse()
    # Calculer le polynôme à factoriser (i.e. G(x) + mu(x))
    A.<x> = PolynomialRing (GF(p))
    Q = A(list (Minv * V(g ** (e - h * d)))) + A(mu)
    # Récupérer les opposés des racines
    beta = [p - Q.roots()[i][0] for i in range (h)]
    # Recouvrer le message clair
    m = [0 for i in range (p)]
    for k in beta :
        m[piInv [alpha.index(k)]] = 1
    return m

def is_not_consistent (pi) :
    size = len(pi) - pi.count(-1)
    s = set(pi) - {-1}
    return size != len(s)

def blapblap (c, p, h, alpha, gpr, r, data) :
    if r ** 2 < h :
        print "r too small, expected at least " + str(int(sqrt(h) + 1))
        return None
    n = h / r
    K = gpr.parent()
    V = K.vector_space()
    pi = [ -1 ] + data + [ -1 for i in range (p - 3) ]
    found = False
    G = [ V(gpr**c[i] - gpr**c[0]) for i in range (n + 1) ]
    M = Matrix(G).transpose()
    X = [ M.solve_right(V(gpr**c[i] - gpr**c[0])) for i in range (p) ]
    u = 0
    while not found and u < p :
        print "u=" + str(u)
        ok = True
        pi = [ -1 ] + data + [ -1 for i in range (p - 3) ]
        i = 3
        while ok and i < p :
            a1 = X[i][1]
            a2 = X[i][2]
            if gcd (int(a2 - u * a1), p) == 1 and a1 != 0 :
                aa = mod (a2 * alpha[pi[2]] - u * a1 * alpha[pi[1]], p) *  mod (a2 - u * a1, p) ** (-1)
                pi[i] = alpha.index(aa)
                if is_not_consistent (pi) :
                    ok = False
            i = i + 1
            print "\t" + str(pi)
        # if ok :
        #     found = True
        u = u + 1
    if not found :
        return "fail"
    else :
        return "MUY BIEN!"

def is_generator (z) :
    K = z.parent()
    q = K.cardinality() - 1
    if z.is_zero() :
        return False
    fac = list(factor(q))
    for div in fac :
        if (z ** (q / div[0])).is_one() :
            return False
    return True

def blyp (PubKey, r) :
    [c, p, h, Q, alpha] = PubKey
    K.<b> = FiniteField(p ** h)
    g = K.random_element()
    while not is_generator(g) :
        g = K.random_element()
    z = g ** ((p ** h - 1) / (p ** r - 1))
    L = []
    for i in range (p ** r - 1) :
        print i
        if gcd (i, p ** r - 1) == 1 :
            L.append(z ** i)
    return L
