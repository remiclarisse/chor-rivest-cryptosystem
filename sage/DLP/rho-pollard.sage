def rho_pollard (h, g) :
    x1, a1, b1 = walk (K.one(), 0, 0, g, h)
    x2, a2, b2 = walk (x1, a1, b1, g, h)
    while x1 != x2 :
        x1, a1, b1 = walk (x1, a1, b1, g, h)
        x2, a2, b2 = walk (x2, a2, b2, g, h)
        x2, a2, b2 = walk (x2, a2, b2, g, h)
    modulo = g.multiplicative_order()
    d = gcd(b1 - b2, modulo)
    l = int(mod ((a2 - a1) / d, modulo / d) * mod ((b1 - b2) / d, modulo / d) ** (-1))
    for i in range (d) :
        l += i * modulo / d
        if g ** l == h :
            return l

def walk (x, a, b, g, h) :
    V = g.parent().vector_space()
    belong = int(sum( V(x).change_ring(IntegerModRing(3)) ))
    modulo = g.multiplicative_order()
    if belong == 1 :
        return g * x, (a + 1) % modulo, b
    elif belong == 0 :
        return x * x, (2 * a) % modulo, (2 * b) % modulo
    else :
        return h * x, a, (b + 1) % modulo
