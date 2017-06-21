def blop (t, p, h) :
    K = t.parent()
    for j in range (h) :
        for a in range (1, p) :
            for b in range (p) :
                if t**(p**j) + a*t + b == 0:
                    print j, a, b
    return None

def blip (c, p, h, alpha, t) :
    r = int(p ** h - 1)
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [ int(mod (log (t + alpha[i], gg), r)) for i in range (p) ]
    C = [ int(mod (c[i]-c[0], r)) for i in range(p) ]
    A = [ int(mod (a[0]-a[i], r)) for i in range(p) ]
    for i in range (p) :
        d = gcd (C[i], r)
        for j in range (p) :
            if A[j] % d == 0 and C[i] != 0 and A[j] != 0:
                gamma = Integer(C[i] / d)
                modulus = Integer(r / d)
                beta = Integer(A[j] / d)
                L = ((gamma).inverse_mod(modulus) * (beta % modulus)) % modulus
                # print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d)
                for w in range (d) :
                    # if w % 100 == 0 and w != 0 :
                    #     print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d) + "\tw:" + str(w)
                    # if ok : sigma(0) = j and sigma(i) = 0
                    D = [ (L * C[l]) % r for l in range(p) ]
                    E = [ (a[k] - a[j]) % r for k in range(p) ]
                    if set(D) == set(E) :
                        return gg**L, int(mod(c[0] - log(t + alpha[j], gg**L), r - 1))
                    L = L + modulus
    return "fail"

def run (p, h) :
    [PubKey, PrivKey] = CRGenerateKeys (p, h);
    [c, p, h, Q, alpha] = PubKey;
    [t, g, sInv, d] = PrivKey;
    gg, dd = blip (c, p, h, alpha, t)
    print "g:" + str(g == gg) + ", d:" + str(d == dd)
