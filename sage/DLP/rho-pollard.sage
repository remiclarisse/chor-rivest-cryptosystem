def rho_pollard (h,g) :
    x1, a1, b1 = walk (h * g, 1, 1, g, h)
    x2, a2, b2 = walk (x1, a1, b1, g, h)
    while x1 != x2 :
        x1, a1, b1 = walk (x1, a1, b1, g, h)
        x2, a2, b2 = walk (x2, a2, b2, g, h)
        x2, a2, b2 = walk (x2, a2, b2, g, h)
    modulo = g.multiplicative_order()
    if gcd(b1 - b2, modulo) == 1 :
        return mod (a2 - a1, modulo) * mod (b1 - b2, modulo) ** (-1)
    else :
        return "Echec"

def walk (x, a, b, g, h) :
    V = g.parent().vector_space()
    belong = int(sum( V(x).change_ring(IntegerModRing(3)) ))
    modulo = g.multiplicative_order()
    if belong == 0 :
        return g * x, (a + 1) % modulo, b
    elif belong == 1 :
        return x * x, (2 * a) % modulo, (2 * b) % modulo
    else :
        return h * x, a, (b + 1) % modulo
