def CRGenerateKeys (q, k) :
    F = FiniteField (q, 'a')
    R = PolynomialRing (F, 'x')
    mu = R.irreducible_element (k)
    K = R.quotient_ring (mu, 'u')
    # Prendre une numérotation de GF(q)
    alpha = F.list()
    shuffle (alpha)
    # Prendre un t algébrique de degré k sur GF(q)
    while True :
        t = K.random_element()
        if t.minpoly().degree() == k :
            break
    # Prendre un g primitif dans GF(q)
    while True :
        g = K.random_element()
        if is_primitive (g) :
            break
    # Prendre une permutation des indices
    s = Permutations([i for i in range (q)]).random_element()
    sInv = [s.index(i) for i in range (q)]
    # Prendre un entier d
    d = randint (0, q ** k - 2)
    # Construire la clé publique
    c = computePubKey (d, t, alpha, s, g, q, k)
    # Retourner les clés !
    PubKey = [c, q, k, alpha]
    PrivKey = [t, g, sInv, d]
    return [PubKey, PrivKey]

def CREncrypt (m, PubKey) :
    # Récuperer la clé publique
    [c, p, h, Q, alpha] = PubKey
    # Retourner le sac à dos
    return mod (sum ([int(m[i])*c[i] for i in range (p)]), p ** h - 1)

def CRDecrypt (e, PubKey, PrivKey) :
    # Récupérer la clé publique
    [c, p, h, Q, alpha] = PubKey
    # Récupérer la clé privée
    [t, g, sInv, d] = PrivKey
    mu = t.minimal_polynomial()
    # Faire le changement de base
    V = t.parent().vector_space()
    Minv = Matrix (GF(p), [V(t ** i) for i in range (h)]).transpose().inverse()
    # Calculer le polynôme à factoriser (i.e. G(x) + mu(x))
    A.<x> = PolynomialRing (GF(p))
    poly = A(list (Minv * V(g ** (e - h * d)))) + A(mu)
    # Récupérer les opposés des racines
    beta = [p - poly.roots()[i][0] for i in range (h)]
    # Recouvrer le message clair
    m = [0 for i in range (p)]
    for k in beta :
        m[sInv [alpha.index(k)]] = 1
    return m

@parallel
def computeLog (h, g, i) :
    return i, log (h, g)

def computePubKey (d, t, alpha, s, g, p, h) :
    w = [ t + alpha[s[i]] for i in range (p) ]
    c = list(computeLog ([ (w[i], g, i) for i in range(p) ]))
    c = [ c[i][1] for i in range (p) ]
    c.sort()
    modulus = p ** h - 1
    c = [ mod (d + c[i][1], modulus) for i in range (p) ]
    return c

def generateRandomMessage (q, k) :
    m = [1 for i in range (k)] + [0 for i in range (q - k)]
    shuffle (m)
    return m

def is_primitive (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise Exception("ring not finite")
    card = card - 1
    for p in prime_factors(card) :
        if x ** (card / p) == A.one() :
            return False
    return True
