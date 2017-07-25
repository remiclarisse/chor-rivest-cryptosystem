def pohlig_hellman (g, h, Fa) :
    n = 1
    for i in range (len (Fa)) :
        n = n * Fa[i][0] ** Fa[i][1]
    K = GF(n + 1)
    a = [0 for i in range (len (Fa))]
    for i in range (len (Fa)) :
        for j in range (1, Fa[i][1] + 1) :
            g0 = g ** (n / (Fa[i][0] ** j))
            h0 = (g0 ** (-a[i])) * (h ** (n / (Fa[i][0] ** j)))
            if h0 != 1 :
                g0 = g ** (n / Fa[i][0])
                b = 1
                t = g0
                while h0 != t :
                    b = b + 1
                    t = t * g0
                a[i] = a[i] + b * Fa[i][0] ** (j - 1)
    moduli = [Fa[i][0] ** Fa[i][1] for i in range (len (Fa))]
    res = CRT (a, moduli)
    return res
