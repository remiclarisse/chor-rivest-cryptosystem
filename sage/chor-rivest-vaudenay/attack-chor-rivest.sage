# reset()
# attach("/home/grace/rclariss/chor-rivest-cryptosystem/sage/chor-rivest.sage")
# attach("/home/grace/rclariss/chor-rivest-cryptosystem/sage/attack-chor-rivest.sage")
# p = 31
# h = 12
# [PubKey, PrivKey] = CRGenerateKeys (p, h)
# m = generateRandomMessage (p, h)
# e = CREncrypt (m, PubKey)
# %time crackedMessage, EquivPrivKey = VaudenayAttack (e, PubKey)

def in_spaned_subspace (v, M) :
    value = True
    try :
        M.solve_right(v)
    except Exception :
        value = False
    finally :
        return value

def firstVaudenayAttack (z, c, p, h, r) :
    n = h / r
    order = p ** r - 1
    K = z.parent()
    V = K.vector_space()
    pourmillage = -1
    maxim = p ** r - 1
    zeta = 1
    i = 1
    while i < p ** r - 1 :
        zeta *= z
        if int(i * (10 ** 3) / maxim) != pourmillage :
            pourmillage = int(i * (10 ** 3) / maxim)
            print "Picking a good generator of GF(" + str(p) + "^" + str(r) + "): \t" + str(pourmillage) + "‰"
        if zeta.multiplicative_order() ==  order :
            G = [ V(zeta**c[w] - zeta**c[0]) for w in range (1, n + 1) ]
            M = Matrix(G).transpose()
            j = 1
            ok = True
            while ok and j < p :
                if not in_spaned_subspace (V(zeta**c[j] - zeta**c[0]), M) :
                    ok = False
                j = j + 1
            if ok :
                print "Picking a good generator of GF(" + str(p) + "^" + str(r) + "): \tdone"
                return zeta
        i += 1
    raise Exception('Unable to find a primitive element gpr in GF(p^r) such that all gpr^ci stands on the same affine sub-space !')

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

def secondVaudenayAttack (c, p, h, alpha, gpr, r, data) : # data === array of size two (e.g. [0, 1])
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

def thirdVaudenayAttack (c, p, h, alpha, pi, gpr, r) :
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

def simplifiedGoldreichAttack (c, p, h, alpha, t, pi) :
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    print "Computing set of logs: ongoing"
    a = [ int(mod (log (t + alpha[pi[i]], gg), r)) for i in range (p) ]
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

# def VaudenayAttack (e, PubKey) :
#     # Récuperer la clé publique
#     [c, p, h, Q, alpha] = PubKey
#     # Choisir un r adéquat
#     for r in divisors (h) :
#         if r ** 2 >= h :
#             break
#     print "Chosen divisor of h: r=" + str(r)
#     K.<b> = FiniteField(p ** h)
#     # Trouver un élément primitif gpr de GF(p^r) tel que les gpr^ci soient dans
#     # le même sous-espace affine
#     g = K.multiplicative_generator()
#     z = g ** ((p ** h - 1) / (p ** r - 1))
#     gpr = firstVaudenayAttack (z, c, p, h, r)
#     # Trouver une permutation d'une clé équivalente
#     alpha = [ K(alpha[i]) for i in range (p) ]
#     pi = secondVaudenayAttack (c, p, h, alpha, gpr, r, [0, 1])
#     piInv = [pi.index(i) for i in range (p)]
#     # Trouver un élément t algébrique de degré h dans GF(p^h)
#     t = thirdVaudenayAttack (c, p, h, alpha, pi, gpr, r)
#     mu = t.minimal_polynomial()
#     # Enfin faire l'attaque de Goldreich, version simplifiée !
#     g, d = simplifiedGoldreichAttack (c, p, h, alpha, t, pi)
#     # Une clé équivalente est construite : on peut déchiffrer le message comme
#     # le destinataire légitime !
#     # Faire le changement de base
#     print "Decrypting the ciphertext"
#     V = t.parent().vector_space()
#     Minv = Matrix (GF(p), [V(t ** i) for i in range (h)]).transpose().inverse()
#     # Calculer le polynôme à factoriser (i.e. G(x) + mu(x))
#     A.<x> = PolynomialRing (GF(p))
#     poly = A(list (Minv * V(g ** (e - h * d)))) + A(mu)
#     # Récupérer les opposés des racines
#     beta = [p - poly.roots()[i][0] for i in range (h)]
#     # Recouvrer le message clair
#     m = [0 for i in range (p)]
#     for k in beta :
#         m[piInv [alpha.index(k)]] = 1
#     print "JOB DONE!!"
#     return m, [t, g, piInv, d]
