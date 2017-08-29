def pohlig_hellman (h, g, fa) :
    n = 1
    for p, i in fa :
        n = n * p ** i
    a = [0 for i in range (len (fa))]
    index = 0
    for p, i in fa :
        for j in range (1, i + 1) :
            g0 = g ** (n / (p ** j))
            h0 = (g0 ** (-a[index])) * (h ** (n / (p ** j)))
            if h0 != 1 :
                g0 = g ** (n / p)
                lg = baby_step_giant_step (g0, h0, p ** j)
                a[index] = a[index] + lg * p ** (j - 1)
        index += 1
    moduli = [ p ** i for p, i in fa ]
    res = CRT (a, moduli)
    return res

def baby_step_giant_step (h, g, n) :
    m = int(ceil (sqrt (n)))
    L = [ g ** i for i in range (m + 1) ]
    u = g ** (-m)
    y = h
    j = 0
    while y not in L :
        y = y * u
        j = j + 1
    return L.index(y) + m * j
