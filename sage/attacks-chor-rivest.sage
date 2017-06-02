def knownGandD (e, PubKey, g, d) :
    c = PubKey[0]
    p = PubKey[1]
    h = PubKey[2]
    alpha = PubKey[4]
    b = [c[i] - d for i in range (p)]
    t = g ** b[0]
    lookUp = [t + alpha[i] for i in range (p)]
    sig = []
    for i in range (p) :
        sig.append(lookUp.index(g ** b[i]))
    sig = [sig.index(i) for i in range (p)]
    # Déchiffrement ordinaire
    mu = t.minimal_polynomial()
    V = t.parent().vector_space()
    Minv = Matrix (GF(p), [V(t ** i) for i in range (h)]).transpose().inverse()
    A.<x> = PolynomialRing (GF(p))
    Q = A(list (Minv * V(g ** (e - h * d)))) + A(mu)
    beta = [p - Q.roots()[i][0] for i in range (h)]
    m = [0 for i in range (p)]
    for k in beta :
        m[sig [alpha.index(k)]] = 1
    return m

def goldreichAttack (e, PubKey, t) :
    c = PubKey[0]
    p = PubKey[1]
    h = PubKey[2]
    alpha = PubKey[4]
    q = p ** h
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [log (t + alpha[i], gg) for i in range (p)]
    A = [mod (a[i] - a[0], q - 1) for i in range (p)]
    C = [mod (c[i] - c[0], q - 1) for i in range (p)]
    i = 0
    while gcd (A[i], q - 1) != 1 and i < p :
        i = i + 1
    if i == p :
        return "Echec : i"
    else :
        aa = mod (A[i], q - 1) ** (-1)
        found = False
        j = 0
        while not found and j < p :
            L = mod (aa * C[j], q - 1)
            if gcd (L, q - 1) == 1 :
                found = True
                for k in range (p) :
                    if L * A[k] not in C :
                        found = False
                        break
                print gg**L
            j = j + 1
    if j == p :
        return "Echec : j", L, gg**L
    return L
