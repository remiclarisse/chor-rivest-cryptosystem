def rho_pollard (p, g, h) :
    K = FiniteField (p)
    global base
    base = g
    global target
    target = h
    global carac
    carac = p
    x1 = 1
    a1 = 0
    b1 = 0
    i = 0
    x1, a1, b1 = walk ([x1, a1, b1])
    x2, a2, b2 = walk ([x1, a1, b1])
    print "-------step "+str(i)+"-------"
    print "x1="+str(x1)+", a1="+str(a1)+", b1="+str(b1)+", i.e "+str(x1)+" = (g**"+str(a1)+") * (h**"+str(b1)+")"
    print "x2="+str(x2)+", a2="+str(a2)+", b2="+str(b2)+", i.e "+str(x2)+" = (g**"+str(a2)+") * (h**"+str(b2)+")"
    while x1 != x2 :
        i = i + 1
        print "-------step "+str(i)+"-------"
        x1, a1, b1 = walk ([x1, a1, b1])
        x2, a2, b2 = walk (walk ([x2, a2, b2]))
        print "x1="+str(x1)+", a1="+str(a1)+", b1="+str(b1)+", i.e "+str(x1)+" = (g**"+str(a1)+") * (h**"+str(b1)+")"
        print "x2="+str(x2)+", a2="+str(a2)+", b2="+str(b2)+", i.e "+str(x2)+" = (g**"+str(a2)+") * (h**"+str(b2)+")"
    print "-------result-------"
    if gcd(b1 - b2, p - 1) == 1 :
        expo = mod (a2 - a1, p - 1) * mod (b1 - b2, p - 1) ** (-1)
        print "exponent="+str(expo)+", i.e "+str(h)+" = "+str(g)+"**"+str(expo)+" (mod "+str(p)+")"
        #return mod (a2 - a1, p - 1) * mod (b1 - b2, p - 1) ** (-1)
    else :
        print "Echec"
        #return "Echec"

def walk (l) :
    x, a, b = l
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
