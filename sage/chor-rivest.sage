### IMPLEMENTATION DU CRYPTOSYSTEME DE CHOR-RIVEST ###

# Les paramètres p et h sont à choisir de sorte que p soit un nombre premier
# et h friable pour pouvoir faire aisement du logarithme discret dans GF(p^h)
# (essayer p = 197 et h = 6)

# Pour générer les clés :
# Keys = CRGenerateKeys (p, h)
# Keys[0] est la clé publique à diffuser et Keys[1] est la clé privée à conserver
# cela peut prendre du temps

# Pour générer aléatoirement un message :
# m = generateRandomMessage (p, h)

# Pour chiffrer un message :
# e = CREncrypt (m, PubKey)
# où PubKey est la clé publique (i.e. Keys[0])

# Pour comparer le déchiffré d'un message à son clair original :
# m == CRDecrypt (e, Keys)
# affiche True si tout s'est bien passé !


def CRGenerateKeys (p, h) :
    q = p ** h
    K.<a> = FiniteField(q)
    Q = a.minimal_polynomial()
    # Prendre une numérotation de GF(p)
    alpha = [K(i) for i in range (p)]
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
    c = [mod (d + log (t + alpha[s[i]], g), q - 1) for i in range (p)]
    # Retourner les clés !
    PubKey = [c, p, h, Q, alpha]
    PrivKey = [t, g, sInv, d]
    return [PubKey, PrivKey]

def CREncrypt (m, PubKey) :
    # Récuperer la clé publique
    c = PubKey[0]
    p = PubKey[1]
    h = PubKey[2]
    q = p ** h
    # Retourner le sac à dos
    return mod (sum ([m[i]*c[i] for i in range (p)]), q-1)

def CRDecrypt (e, Keys) :
    # Récupérer la clé publique
    p = Keys[0][1]
    h = Keys[0][2]
    alpha = Keys[0][4]
    # Récupérer la clé privée
    t = Keys[1][0]
    mu = t.minimal_polynomial()
    g = Keys[1][1]
    sInv = Keys[1][2]
    d = Keys[1][3]
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
        m[sInv [alpha.index(k)]] = 1
    return m

def generateRandomMessage (p, h) :
    m = [1 for i in range (h)] + [0 for i in range (p - h)]
    shuffle (m)
    return m

def transformMessage (m, p, h) :
    b = floor( log (binomial (p, h),2))
    n = []
    while len (m) > 0 :
        k = m[:b]
        m = m[b:]
        n.append(int(k,2))
    for i in range (len (n)) :
        y = ""
        k = h
        for j in range (p) :
            if n[i] >= binomial (p - (j + 1), k) :
                y = y + '1'
                n[i] = n[i] - binomial (p - (j + 1), k)
                k = k - 1
            else :
                y = y + '0'
        n[i] = y
    return n

def recoverMessage (e, p, h) :
    a = []
    for i in range(len (e)) :
        n = 0
        k = h
        for j in range (p) :
            if e[i][j] == '1' :
                n = n + binomial(p - (j + 1), k)
                k = k - 1
        a.append(n)
    m = ""
    b = floor( log (binomial (p, h),2))
    for n in a :
        s = bin (n)[2:]
        l = len (s)
        if n != a[-1] :
            for i in range (b - l) :
                s = '0' + s
        m = m + s
    return m

# encodeMessage("Why did the topologist's marriage fail ? Because he thought that arbitrary unions were open !")
def encodeMessage (m) :
    import binascii
    hexa = binascii.hexlify(m)
    bina = ""
    for chara in hexa :
        bina = bina + '{0:08b}'.format(Integer(chara,16))
    return bina

def decodeMessage (c) :
    import binascii
    hexa = ""
    while len (c) > 0 :
        bina = c[:8]
        c = c[8:]
        hexa = hexa + '{0:01x}'.format(Integer(bina,2))
    return binascii.unhexlify(hexa)
