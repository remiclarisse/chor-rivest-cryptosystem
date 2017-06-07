def goldreichAttack (e, PubKey, t) :
    c = PubKey[0]
    p = PubKey[1]
    h = PubKey[2]
    alpha = PubKey[4]
    q = p ** h
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [log (t + alpha[i], gg) for i in range (p)]
    found = False
    i0 = 0
    while not found and i0 < p:
        j0 = i0 + 1
        while not found and j0 < p:
            if gcd (c[j0] - c[i0], p ** h - 1) == 1 :
                found = True
            else :
                j0 = j0 + 1
    if not found :
        return "ECHEC : 1"
    found = False
    b = 0
    while not found and b < p : # suppose sigma(i0) = b
        for i in range (p) :
            l = mod (log (t + alpha[i], gg) - log (t + alpha[b] ,gg), p ** h - 1) * mod (c[j0] - c[i0], p ** h - 1) ** (-1)
            qsd = [mod (log (t + alpha[j], gg) - log (t + alpha[b], gg), p ** h - 1) for j in range (p)]
            wxc = [mod (l * (c[j] - c[i0]), p ** h - 1) for j in range (p)]
            if set(qsd) == set(wxc) :
                found = True
                return gg ** l
        b = b + 1
    return "ECHEC : 2"
