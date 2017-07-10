# p = 2^15 * 3^6 * 5^4 + 1 ;
# g = GF(p)(4107883536) ;
# h = GF(p)(572067100) ;
# Fa=[[2,15],[3,6],[5,4]] ;
# log(h,g) = 2074276358 ;

def PohligHellmanField (g, h, Fa) :
    n = 1
    for i in range (len (Fa)) :
        n = n * Fa[i][0] ** Fa[i][1]
    K = GF(n + 1)
    G = [[K(g ** (n / (Fa[i][0] ** j))) for j in range (1, Fa[i][1] + 1)] for i in range (len (Fa))]
    H = [[K(h ** (n / (Fa[i][0] ** j))) for j in range (1, Fa[i][1] + 1)] for i in range (len (Fa))]
    a = [0 for i in range (len (Fa))]
    for i in range (len (Fa)) :
        for j in range (Fa[i][1]) :
            g0 = G[i][j]
            h0 = H[i][j]
            u = g0 ** (-a[i])
            h0 = h0 * u
            if h0 != 1 :
                g0 = G[i][0]
                b = 1
                T = g0
                while h0 != T :
                    b = b + 1
                    T = T * g0
                a[i] = a[i] + b * Fa[i][0] ** j
    moduli = [Fa[i][0] ** Fa[i][1] for i in range (len (Fa))]
    res = CRT (a, moduli)
    return res

def PohligHellmanFieldBis (g, h, Fa) :
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
