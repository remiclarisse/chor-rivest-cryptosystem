def blop (t, p, h) :
    K = t.parent()
    for j in range (h) :
        for a in range (1, p) :
            for b in range (p) :
                if t**(p**j) + a*t + b == 0:
                    print j, a, b
    return None

def blip (c, p, h, alpha, t) :
    r = p ** h - 1
    K = t.parent()
    gg = K.multiplicative_generator()
    a = [ mod(log (t + alpha[i], gg), r) for i in range (p) ]
    C = [ mod (c[i]-c[0], r) for i in range(p) ]
    A = [ mod (a[0]-a[i], r) for i in range(p) ]
    for i in range (p) :
        if gcd (C[i], r) == 1 :
            for j in range (p) :
                L = mod (C[i], r) ** (-1) * mod (A[j], r)
                # if ok : sigma(0) = j and sigma(i) = 0
                still_looking = True
                k = 0
                D = [ mod (L * C[l], r) for l in range(p) ]
                while still_looking and k < p:
                    if mod (a[k] - a[j], r) not in D :
                        still_looking = False
                    else :
                        k = k + 1
                if k == p :
                    return gg**L, mod(c[0] - log(t + alpha[j], gg**L), r-1)
        else :
            for j in range (p) :
                if gcd (C[i], r) != 0 :
                    if A[j] % gcd (C[i], r) == 0 :
                        d = gcd (C[i], r)
                        for w in range (d) :
                            L = mod (C[i] / d, r / d) ** (-1) * mod (A[j] / d, r / d) + w * (r / d)
                            # if ok : sigma(0) = j and sigma(i) = 0
                            still_looking = True
                            k = 0
                            D = [ mod (L * C[l], r) for l in range(p) ]
                            while still_looking and k < p:
                                if mod (a[k] - a[j], r) not in D :
                                    still_looking = False
                                else :
                                    k = k + 1
                            if k == p :
                                return gg**L, mod(c[0] - log(t + alpha[j], gg**L), r-1)
    return "fail"
