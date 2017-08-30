def pohlig_hellman (h, g, fa) :
    n = 1
    for p, i in fa :
        n = n * p ** i
    a = [0 for i in range (len (fa))]
    index = 0
    for p, i in fa :
        g0 = g ** (n / p)
        for j in range (1, i + 1) :
            h0 = (  g ** ((n  / (p ** j)) * (-a[index])) * h ** (n / (p ** j))  )
            if h0 != h0.parent().one() :
                lg = baby_step_giant_step (g0, h0, p)
                a[index] = a[index] + lg * p ** (j - 1)
        index += 1
    moduli = [ p ** i for p, i in fa ]
    res = CRT (a, moduli)
    return res

def baby_step_giant_step (h, g, n) :
    m = int(ceil (sqrt (n)))
    L = []
    for i in range (m + 1) :
        u = g ** i
        L += [ hash(u) ]
    u = u ** (-1)
    y = h
    j = 0
    while hash(y) not in L :
        y = y * u
        j = j + 1
    return L.index(hash(y)) + m * j
