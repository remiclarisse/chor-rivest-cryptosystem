### IMPLEMENTATION DU CRYPTOSYSTEME DE CHOR-RIVEST ###

# Les paramètres p et h sont à choisir de sorte que p soit un nombre premier
# et h friable pour pouvoir faire aisement du logarithme discret dans GF(p^h)
# (essayer p = 197 et h = 24)

# Pour générer les clés :
# [PubKey, PrivKey] = CRGenerateKeys (p, h)
# PubKey est la clé publique à diffuser et PrivKey est la clé privée à conserver
# L'exécution de cette commande peut prendre du temps

# Pour générer aléatoirement un message :
# m = generateRandomMessage (p, h)

# Pour chiffrer un message :
# e = CREncrypt (m, PubKey)

# Pour comparer le déchiffré d'un message à son clair original :
# m == CRDecrypt (e, PubKey, PrivKey)
# affiche True si tout s'est bien passé !

def CRGenerateKeys (p, h) :
    q = p ** h
    K.<a> = FiniteField(q)
    Q = a.minimal_polynomial()
    # Prendre une numérotation de GF(p)
    alpha = [int(i) for i in range (p)]
    shuffle (alpha)
    # Prendre un t algébrique de degré h sur GF(p)
    while True :
        t = K.random_element()
        if t.minimal_polynomial().degree() == h :
            break
    # Prendre un g primitif dans GF(q)
    while True :
        g = K.random_element()
        if g.multiplicative_order() == q - 1 :
            break
    # Prendre une permutation des indices
    s = Permutations([i for i in range (p)]).random_element()
    sInv = [s.index(i) for i in range (p)]
    # Prendre un entier d
    d = randint (0, q - 2)
    # Construire la clé publique
    c = computePubKey (d, t, alpha, s, g, p, h)
    # Retourner les clés !
    PubKey = [c, p, h, Q, alpha]
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
    return i, discrete_log (h, g)

def computePubKey (d, t, alpha, s, g, p, h) :
    w = [ t + alpha[s[i]] for i in range (p) ]
    c = list(computeLog ([ (w[i], g, i) for i in range(p) ]))
    c = [ c[i][1] for i in range (p) ]
    c.sort()
    modulus = p ** h - 1
    c = [ mod (d + c[i][1], modulus) for i in range (p) ]
    return c

def generateRandomMessage (p, h) :
    m = [1 for i in range (h)] + [0 for i in range (p - h)]
    shuffle (m)
    return m
