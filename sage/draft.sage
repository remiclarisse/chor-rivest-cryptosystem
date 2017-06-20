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
        print "i:" + str(i)
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
                        print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d)
                        for w in range (d) :
                            if w % 100 == 0 and w != 0 :
                                print "i:" + str(i) + "\tj:" + str(j) + "\td:" + str(d) + "\tw:" + str(w)
                            L = mod (int(C[i]) / int(d), int(r) / int(d)) ** (-1) * mod (int(A[j]) / int(d), int(r) / int(d)) + w * (int(r) / int(d))
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
