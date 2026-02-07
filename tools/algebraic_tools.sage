def construct_quadratic_system(R, n, m, planted=False):
    fs = [sum(x * R.base_ring().random_element() for x in R.monomials_of_degree(2)) for _ in range(m)]
    sol = [0 for _ in range(n)]
    if planted:
        # We assume x_n != 0
        x_n = R.gens()[-1]
        sol = [R.base_ring().random_element() for _ in range(n-1)] + [1]
        fs = [f - f(sol)*x_n*x_n for f in fs]
        assert [f(sol) == 0 for f in fs]
    return fs, sol

def construct_Macaulay_matrix(R, d, fs):
    m = len(fs)

    mons_mac = R.monomials_of_degree(d)[::-1]
    mults_mac = R.monomials_of_degree(d-2)[::-1]
    mac_entries = {(m * j + i, mons_mac.index(mon)): coeff for j, mult in enumerate(mults_mac) for i, f in enumerate(fs) for (coeff, mon) in mult * f if mon in mons_mac}
    return matrix(R.base_ring(), m * len(mults_mac), len(mons_mac), mac_entries)

def construct_relation_matrix(R, d, fs):
    m = len(fs)

    mons_rel = R.monomials_of_degree(d-2)[::-1]
    mults_rel = R.monomials_of_degree(d-4)[::-1]
    rel_entries = {}
    block = 0
    for i in range(m):
        for j in range(i+1, m):
            rel_entries.update({(u * binomial(m, 2) + block, m * mons_rel.index(mon) + i): coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[j]})
            rel_entries.update({(u * binomial(m, 2) + block, m * mons_rel.index(mon) + j): -coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[i]})
            block += 1
    return matrix(R.base_ring(), binomial(m, 2) * len(mults_rel), m * len(mons_rel), rel_entries)
