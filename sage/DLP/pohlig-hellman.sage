def pohlig_hellman (h, g) :
    n = g.multiplicative_order()
    fa = list(factor (n))
    a = [0 for i in range (len (fa))]
    for p, i in fa :
        g0 = g ** (n / p)
        ind = fa.index((p, i))
        for j in range (1, i + 1) :
            h0 = (  h * g ** (-a[ind])  ) ** (n / (p ** j))
            b = baby_step_giant_step (h0, g0, p)
            a[ind] = a[ind] + b * p ** (j - 1)
    moduli = [ p ** i for p, i in fa ]
    res = CRT (a, moduli)
    return res

def baby_step_giant_step (h, g, n) :
    m = int(ceil (sqrt (n)))
    L = []
    for i in range (m + 1) :
        u = g ** i
        L += [ u ]
    u = u ** (-1)
    y = h
    j = 0
    while y not in L :
        y = y * u
        j = j + 1
    return L.index(y) + m * j
