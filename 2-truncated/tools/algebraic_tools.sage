def a_sequences(n, d):
    if d == 0: yield tuple()
    for i in range(n, 0, -1):
        for a in a_sequences(i-1, d-1):
            yield a + (i,)

def bitruncated_monomials(R,d):
    if d == 0: return [R.one()]
    return [product([R.gen(i-1) for i in indices]) for indices in a_sequences(len(R.gens()),d)]  

def construct_quadratic_system(R, n, m, o=0):
    v = n-o
    fs = []
    while True:
        transform = random_matrix(R.base_ring(),n,n)
        if det(transform) != 0: break
    for _ in range(m):
        centralMap = block_matrix([[random_matrix(R.base_ring(),v,v),random_matrix(R.base_ring(),v,o)],[random_matrix(R.base_ring(),o,v),zero_matrix(R.base_ring(),o,o)]])
        F = transform.transpose() * centralMap * transform
        fs.append( sum([R.gen(i) * R.gen(j) * F[i,j] for i in range(n) for j in range(n) if i != j]) )
    return fs, transform

def construct_Macaulay_matrix(R, d, fs):
    m = len(fs)

    mons_mac = bitruncated_monomials(R,d)[::-1]
    mults_mac = bitruncated_monomials(R,d-2)[::-1]
    mac_entries = {(m * j + i, mons_mac.index(mon)): coeff for j, mult in enumerate(mults_mac) for i, f in enumerate(fs) for (coeff, mon) in mult * f if mon in mons_mac}
    return matrix(R.base_ring(), m * len(mults_mac), len(mons_mac), mac_entries)

def construct_relation_matrix(R, d, fs):
    m = len(fs)

    mons_rel = bitruncated_monomials(R,d-2)[::-1]
    mults_rel = bitruncated_monomials(R,d-4)[::-1]
    rel_entries = {}
    block = 0
    for i in range(m):
        for j in range(i, m):
            rel_entries.update({(u * binomial(m, 2) + block, m * mons_rel.index(mon) + i): coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[j] if mon in mons_rel})
            rel_entries.update({(u * binomial(m, 2) + block, m * mons_rel.index(mon) + j): -coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[i] if mon in mons_rel})
            block += 1
    return matrix(R.base_ring(), binomial(m+1, 2) * len(mults_rel), m * len(mons_rel), rel_entries)