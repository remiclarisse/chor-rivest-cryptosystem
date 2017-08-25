def rho_pollard (g, h) :
    K = FiniteField (p)
    x1, a1, b1 = walk (x1, a1, b1)
    x2, a2, b2 = walk (x1, a1, b1)
    while x1 != x2 :
        i = i + 1
        x1, a1, b1 = walk ([x1, a1, b1])
        x2, a2, b2 = walk (walk ([x2, a2, b2]))
    if gcd(b1 - b2, p - 1) == 1 :
        expo = mod (a2 - a1, p - 1) * mod (b1 - b2, p - 1) ** (-1)
        #return mod (a2 - a1, p - 1) * mod (b1 - b2, p - 1) ** (-1)
    else :
        print "Echec"
        #return "Echec"

def walk (x, a, b) :
    if 0 < x and x <= carac / 3 :
        y = (base * x) % (carac)
        c = (a + 1) % (carac - 1)
        d = (b) % (carac - 1)
    elif  carac / 3 < x and x <= 2 * carac / 3 :
        y = (x * x) % (carac)
        c = (2 * a) % (carac - 1)
        d = (2 * b) % (carac - 1)
    else :
        y = (target * x) % (carac)
        c = (a) % (carac - 1)
        d = (b + 1) % (carac - 1)
    return [y, c, d]
