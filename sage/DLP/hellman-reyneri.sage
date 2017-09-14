def hellman_reyneri (g, h) :
    K = g.parent()
    R = K.base()
    F = R.base_ring()
    q = F.cardinality()
    card = K.cardinality()
    k = log (card, q)
    card -= 1
    # collecte de relations
    print "Collecte de relations"
    d = int(ceil (sqrt (0.5 * k * log (k, q))))
    rs, rel_liste = make_relations (g, d, q ** (d + 1), card)
    print "Phase d'algèbre linéaire"
    basis, logs_basis = make_smoothness_basis (rel_liste, rs, card)
    # descente individuelle
    print "Descente individuelle"
    logs = [ mod (descent_log (elem_h, g, basis, logs_basis, card), card) for elem_h in h ]
    return logs

def descent_log (elem_target, gen, basis, logs_basis, card) :
    P = elem_target.lift()
    r = 0
    while not is_in_basis (P, basis) :
        r = randint (1, card - 1)
        P = (elem_target * gen ** r).lift()
    lo = 0
    if P.leading_coefficient() != P.parent().one() :
        lo += logs_basis[basis.index(P.leading_coefficient())]
    for poly, mult in list(P.factor()) :
        lo += mult * logs_basis[basis.index(poly)]
    return lo - r

def is_in_basis (x, basis) :
    for poly, mult in list(x.factor()) :
        if poly not in basis :
            return False
    return True

def make_smoothness_basis (rel_liste, rs, modulus) :
    A, basis = make_matrix_relation_n_get_basis (rel_liste, modulus)
    B = vector (IntegerModRing(modulus), rs)
    return basis, A.solve_right(B)

def make_matrix_relation_n_get_basis (liste, modulus) :
    unknowns = []
    for lc, fa in liste :
        if lc not in unknowns and lc != lc.parent().one() :
            unknowns.append(lc)
        for poly, mult in fa :
            if poly not in unknowns :
                unknowns.append(poly)
    M = []
    for lc, fa in liste :
        lig = [ 0 for i in range (len (unknowns)) ]
        if lc in unknowns :
            lig[unknowns.index(lc)] = 1
        for poly, mult in fa :
            lig[unknowns.index(poly)] = mult
        M = M + [lig]
    M = Matrix(IntegerModRing(modulus), M, sparse=True)
    return M, unknowns

def make_relations (gen, threshold, numb, card) :
    sol = []
    rel_liste = []
    tirage = 0
    while tirage < numb :
     r = randint (1, card - 1)
     P = (gen ** r).lift()
     if P.degree() <= threshold :
         lc = P.leading_coefficient()
         P = P * lc ** (-1)
         tirage = tirage + 1
         if r not in sol :
             sol.append(r)
             rel_liste.append((lc, list(P.factor())))
    return sol, rel_liste

def make_change_basis_matrix (q, k, F, R, u, Ik) :
    i = 0
    f = [ F(0) for i in range (k) ]
    while (R(f)(u)).minpoly() != Ik and i < q**k :
        i += 1
        f = next_poly(f)
    if i == q ** k :
        raise Exception("can't find f such that v=f(u)")
    t = R(f)(u)
    return Matrix( [ list(t**i) for i in range (k) ] ).transpose()

def pick_primitive_element (K) :
    x = K.random_element()
    while not is_primitive (x) :
        x = K.random_element()
    return x

def is_primitive (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise Exception("ring not finite")
    card = card - 1
    for p in prime_factors(card) :
        if x ** (card / p) == A.one() :
            return False
    return True

def mult_order (x) :
    A = x.parent()
    card = A.cardinality()
    if card.parent() != ZZ :
        raise Exception("ring not finite")
    card = card - 1
    for d in divisors(card) :
        if x ** d == A.one() :
            return d

def next_poly (P) :
    F = P[0].parent()
    if P[0] == 1 :
        return [F(0)] + P[1:]
    elif F.next(P[0]) != 1 :
        return [F.next(P[0])] + P[1:]
    else :
        return [F(1)] + next_poly(P[1:])
