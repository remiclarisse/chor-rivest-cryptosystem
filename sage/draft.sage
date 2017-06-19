def blop (t, p, h) :
    K = t.parent()
    for j in range (h) :
        for a in range (1, p) :
            for b in range (p) :
                if t**(p**j) + a*t + b == 0:
                    print j, a, b
    return None
