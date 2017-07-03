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

def blup (c, p, h, alpha, s, gpr, r) :
    n = h / r
    K = gpr.parent()
    A.<X> = PolynomialRing(K)
    Q = A(0)
    for i in range (n + 1) :
        L = A(1)
        for k in range (n + 1) :
            if k != i:
                L = L * A((alpha[s[i]] - alpha[s[k]]) ** (-1) * (X - alpha[s[k]]))
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
        print newnbs # verbose mode
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
