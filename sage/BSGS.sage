# K.<a> = GF(2 ** 19)
# g = K(a)
# h = K(a^17 + a^13 + a^12 + a^9 + a^6 + a^5 + a^3 + 1)
# 248282 = log (h, g)

def BSGS (g, h, r) :
    m = int(ceil (sqrt (r)))
    L = [g ** i for i in range (m + 1)]
    u = g ** (-m)
    y = h
    j = 0
    while y not in L :
        y = y * u
        j = j + 1
    return L.index(y) + m * j
