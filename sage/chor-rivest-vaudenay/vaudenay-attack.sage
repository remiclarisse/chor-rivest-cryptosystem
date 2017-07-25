def VaudenayAttack (PubKey) :
    # Récuperer la clé publique
    [c, p, h, Q, alpha] = PubKey
    # Choisir un r adéquat
    for r in divisors (h) :
        if r ** 2 >= h :
            break
    print "Chosen divisor of h: r=" + str(r) + ", divisors are: " + str(divisors(r))
    # Trouver un élément primitif gpr de GF(p^r) tel que les gpr^ci soient dans
    # le même sous-espace affine
    print "Computing primitive element gpr..."
    gpr = firstVaudenayAttack (PubKey, r)
    K = gpr.parent()
    # Trouver une permutation d'une clé équivalente
    PubKey[4] = [ K(alpha[i]) for i in range (p) ]
    print "Computing permutation pi..."
    pi = secondVaudenayAttack (PubKey, gpr, r, [0, 1])
    piInv = [pi.index(i) for i in range (p)]
    # Trouver un élément t algébrique de degré h dans GF(p^h)
    print "Computing algebraic element t..."
    t = thirdVaudenayAttack (PubKey, pi, gpr, r)
    # Enfin faire l'attaque de Goldreich, version simplifiée !
    print "Computing primitive element g and integer d..."
    g, d = simplifiedGoldreichAttack (PubKey, t, pi)
    print "JOB DONE"
    return [t, g, piInv, d]

def firstVaudenayAttack (PubKey, r) :
    ri = divisors(r)
    K.<b> = FiniteField(p ** h)
    gpri = [ compute_gp (PubKey, K) ] + [ K(1) for i in range (1, len (ri)) ]
    for i in range (1, len (ri)) :
        dag = find_DAG(ri[i])
        gpri[i] = compute_gpri (PubKey, ri[i], [ ri[j] for j in range (len(ri)) if ri[j] in dag ], [ gpri[j] for j in range (len (gpri)) if ri[j] in dag ])
    if gpri[-1] != 1 :
        return gpri[-1]
    else :
        raise Exception('Unable to find a primitive element gpr in GF(p^r) such that all gpr^ci stands on the same affine sub-space !')

def compute_gp (PubKey, K) :
    [c, p, h, Q, alpha] = PubKey
    g = K.multiplicative_generator()
    gamma = g ** ((p ** h - 1) / (p  - 1)) # élément primitif de GF(p)
    i = 0
    while i < p :
        if gcd (i, p - 1) == 1 :
            z = gamma ** i
            e = 1
            ok = True
            while e < (p - 1) / h and ok :
                res = 0
                for j in range (p) :
                    res += z ** (e * c[j])
                if res != 0 :
                    ok = False
                e += 1
            if ok :
                return K(z)
        i += 1
    raise Exception('Unable to find a primitive element gpr in GF(p) such that all gpr^ci stands on the same affine sub-space !')

def compute_gpri (PubKey, r, ri, gpri) :
    [c, p, h, Q, alpha] = PubKey
    n = h / r
    K = gpri[0].parent()
    g = K.multiplicative_generator()
    gamma = g ** ((p ** h - 1) / (p ** r - 1)) # élément primitif de GF(p^r)
    a = []
    for i in range (len (ri)) :
        l = log (gpri[i], gamma)
        d = 0
        for k in range (int(r / ri[i])) :
            d += p ** (k * ri[i])
        a.append((l / d) % (p ** ri[i] - 1))
    modulii = [ p ** ri[i] - 1 for i in range (len (ri)) ]
    resCRT = CRT(a, modulii)
    beta = gamma ** resCRT
    ell = lcm (modulii)
    x = 0
    while x < ((p ** r - 1) / ell) :
        e = 1
        ok = True
        while e < (p - 1) * r / h and ok :
            res = 0
            for i in range (p) :
                res += beta ** (e * c[i]) * gamma ** (e * c[i] * ell * x)
            if res != 0 :
                ok = False
            e += 1
        if ok :
            return beta * gamma ** (ell * x)
        x += 1

def find_DAG (r) :
    if r.is_prime() or r == 1:
        return [1]
    divs = []
    for k, i in list (factor (r))[::-1] :
        if k != 1 and k != r :
            divs.append(r / k)
    return divs

def in_spaned_subspace (v, M) :
    value = True
    try :
        M.solve_right(v)
    except Exception :
        value = False
    finally :
        return value

def secondVaudenayAttack (PubKey, gpr, r, data) : # data === array of size two (e.g. [0, 1])
    [c, p, h, Q, alpha] = PubKey
    if r ** 2 < h :
        raise Exception('r too small')
    n = h / r
    K = gpr.parent()
    V = K.vector_space()
    pi = [ -1 ] + data + [ -1 for i in range (p - 3) ]
    G = [ V(gpr**c[i] - gpr**c[0]) for i in range (n + 1) ]
    M = Matrix(G).transpose()
    X = [ M.solve_right(V(gpr**c[i] - gpr**c[0])) for i in range (p) ]
    u = 0
    while u < p :
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
        if ok :
            I = [ i for i in range (p) if pi[i] == -1 ]
            S = [ i for i in range (p) if i not in pi ]
            count  = len (I)
            G = Permutations([ i for i  in range (count) ])
            for el in G :
                for i in range (count) :
                    pi[I[i]] = S[el[i]]
                i = 0
                sure = True
                while sure and i < p :
                    Q = K(0)
                    for j in range (n + 1) :
                        L = K(1)
                        for k in range (n + 1) :
                            if k != j :
                                L = L * (alpha[pi[j]] - alpha[pi[k]]) ** (-1) * (alpha[pi[i]] - alpha[pi[k]])
                        Q = Q + gpr ** c[j] * L
                    if gpr ** c[i] != Q :
                        sure = False
                    i += 1
                if sure :
                    return pi
        u = u + 1
    raise Exception('Unable to find a permutation!')

def is_not_consistent (pi) :
    size = len(pi) - pi.count(-1)
    s = set(pi) - {-1}
    return size != len(s)

def thirdVaudenayAttack (PubKey, pi, gpr, r) :
    [c, p, h, Q, alpha] = PubKey
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
    if len(R) == 0 :
        raise Exception('Unable to find an element t!')
    else :
        return -R[0][0]

def simplifiedGoldreichAttack (PubKey, t, pi) :
    [c, p, h, Q, alpha] = PubKey
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    a = computePubKey (0, t, alpha, pi, gg, p, h)
    C = [ Integer(mod (c[i]-c[0], r)) for i in range(p) ]
    A = [ int(mod (a[i]-a[0], r)) for i in range(p) ]
    for i in range (p) :
        if gcd (C[i], r) == 1 :
            L = (C[i]).inverse_mod(r) * (A[i])
            D = [ (L * C[l]) % r for l in range(p) ]
            if set(D) == set(A) :
                return gg**L, int(mod(c[0] - log(t + alpha[pi[0]], gg**L), r - 1))
    raise Exception('Unable to find primitive element g and integer d!')
