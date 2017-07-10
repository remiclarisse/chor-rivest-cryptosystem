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
    gpr = firstVaudenayAttack (PubKey, r)
    # Trouver une permutation d'une clé équivalente
    alpha = [ K(alpha[i]) for i in range (p) ]
    pi = secondVaudenayAttack (PubKey, gpr, r, [0, 1])
    piInv = [pi.index(i) for i in range (p)]
    # Trouver un élément t algébrique de degré h dans GF(p^h)
    t = thirdVaudenayAttack (PubKey, pi, gpr, r)
    print "Computing minimal polynomial: ongoing"
    mu = t.minimal_polynomial()
    print "Computing minimal polynomial: done"
    # Enfin faire l'attaque de Goldreich, version simplifiée !
    g, d = simplifiedGoldreichAttack (PubKey, t, pi)
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
    pourmillage = -1
    maxim = p  - 1
    i = 0
    while i < p :
        if int(i * (10 ** 3) / maxim) != pourmillage :
            pourmillage = int(i * (10 ** 3) / maxim)
            print "Picking a good generator of GF(" + str(p) + "): \t" + str(pourmillage) + "‰"
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
                print "Picking a good generator of GF(" + str(p) + "): \tdone"
                return K(z)
        i += 1
    return "fail"

def compute_gpri (PubKey, r, ri, gpri) :
    [c, p, h, Q, alpha] = PubKey
    n = h / r
    K = gpri[0].parent()
    g = K.multiplicative_generator()
    gamma = g ** ((p ** h - 1) / (p ** r - 1)) # élément primitif de GF(p^r)
    a = []
    print "Building an almost generator of GF(" + str(p) + "^" + str(r) + ") with CRT: ongoing"
    for i in range (len (ri)) :
        l = log (gpri[i], gamma)
        d = 0
        for k in range (int(r / ri[i])) :
            d += p ** (k * ri[i])
        a.append((l / d) % (p ** ri[i] - 1))
    modulii = [ p ** ri[i] - 1 for i in range (len (ri)) ]
    resCRT = CRT(a, modulii)
    print "Building an almost generator of GF(" + str(p) + "^" + str(r) + ") with CRT: done"
    beta = gamma ** resCRT
    ell = lcm (modulii)
    pourmillage = -1
    maxim = ((p ** r - 1) / ell) - 1
    x = 0
    while x < ((p ** r - 1) / ell) :
        if int(x * (10 ** 3) / maxim) != pourmillage :
            pourmillage = int(x * (10 ** 3) / maxim)
            print "Making a good generator of GF(" + str(p) + "^" + str(r) + "): \t" + str(pourmillage) + "‰"
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
            print "Making a good generator of GF(" + str(p) + "^" + str(r) + "): \tdone"
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

def secondVaudenayAttack (PubKey, gpr, r, data) : # data === array of size two (e.g. [0, 1])
    [c, p, h, Q, alpha] = PubKey
    n = h / r
    K = gpr.parent()
    newnbs = data + [i for i in range (2, n)] + [n - 1]
    found = False
    ok = True
    pourmillage = -1
    maxim = binomial(p - 2, n - 1) * factorial (n - 1)
    tracker = 1
    while not found and ok:
        if int(tracker * 1000 / maxim) != pourmillage :
            pourmillage = int(tracker * 1000 / maxim)
            print "Picking a permutation such that π(0)=" + str(data[0]) + " and π(1)=" + str(data[1]) + " : \t" + str(pourmillage) + "‰"
        ok = increment (newnbs, data[1], p)
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
        tracker += 1
    if not ok :
        raise Exception('Unable to find a permutation!')
    else :
        print "Picking a permutation such that π(0)=" + str(data[0]) + " and π(1)=" + str(data[1]) + " : \tdone"
        return sig

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
    print "Picking an algebraic element of GF(" + str(p) + "^" + str(h) + ") of degree " + str(h)
    if len(R) == 0 :
        raise Exception('Unable to find an element t!')
    else :
        return -R[0][0]

def simplifiedGoldreichAttack (PubKey, t, pi) :
    [c, p, h, Q, alpha] = PubKey
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    print "Computing set of logs: ongoing"
    a = computePubKey (0, t, alpha, pi, gg, p, h)
    print "Computing set of logs: done"
    C = [ Integer(mod (c[i]-c[0], r)) for i in range(p) ]
    A = [ int(mod (a[i]-a[0], r)) for i in range(p) ]
    pourmillage = -1
    maxim = p
    tracker = 1
    for i in range (p) :
        if int(i * 1000 / maxim) != pourmillage :
            pourmillage = int(i * 1000 / maxim)
            print "Finding the generator g of GF(" + str(p) + "^" + str(h) + ") and the integer d: \t" + str(pourmillage) + "‰"
        if gcd (C[i], r) == 1 :
            L = (C[i]).inverse_mod(r) * (A[i])
            D = [ (L * C[l]) % r for l in range(p) ]
            if set(D) == set(A) :
                print "Finding the generator g of GF(" + str(p) + "^" + str(h) + ") and the integer d: \tdone"
                return gg**L, int(mod(c[0] - log(t + alpha[pi[0]], gg**L), r - 1))
        tracker += 1
    raise Exception('Unable to find primitive element g and integer d!')
