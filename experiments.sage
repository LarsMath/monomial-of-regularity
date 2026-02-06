import math
from prediction_tools import *

def construct_quadratic_system(R, n, m, planted):
    fs = [sum(x * R.base_ring().random_element() for x in R.monomials_of_degree(2)) for _ in range(m)]
    if planted:
        # We assume x_n != 0
        xn = R.gens()[-1]
        sol = [R.base_ring().random_element() for _ in range(n-1)] + [1]
        fs = [f - f(sol)*xn*xn for f in fs]
        assert [f(sol) == 0 for f in fs]
    return fs

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
            rel_entries.update({(u * math.comb(m, 2) + block, m * mons_rel.index(mon) + i): coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[j]})
            rel_entries.update({(u * math.comb(m, 2) + block, m * mons_rel.index(mon) + j): -coeff for u, mult in enumerate(mults_rel) for (coeff, mon) in mult * fs[i]})
            block += 1
    return matrix(R.base_ring(), math.comb(m, 2) * len(mults_rel), m * len(mons_rel), rel_entries)


# In this method we verify that random systems are monomial semi-regular
# For systems with a solution we verify that they are regular up to nullity 1
# We do this by finding the lowest mu with a non-trivial kernel vector
# At this point, the nullity should be exactly 0 or 1
def verify_monomial_semi_regularity(n, m, planted, q=251):
    # Note that we dehomogenize!
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order='degrevlex')
    dreg = degree_of_regularity(n, m, solution=planted)

    # Construct system
    fs = construct_quadratic_system(R, n, m, planted)

    # Construct Macaulay matrix
    mac = construct_Macaulay_matrix(R, dreg, fs)

    # Pivots of the actual syzygies
    ker_mac_pivots = mac.transpose().right_kernel_matrix().pivots()

    if dreg >= 4:
        rel = construct_relation_matrix(R, dreg, fs)
        # Pivots of the trivial syzygies
        rel_pivots = rel.pivots()
    else:
        rel_pivots = []

    if all(p in rel_pivots for p in ker_mac_pivots):
        # All syzygies are trivial
        return True
    
    first_non_trivial = next(p for p in ker_mac_pivots[::-1] if p not in rel_pivots)

    syzygies = sum(1 for p in ker_mac_pivots if p >= first_non_trivial)
    rows = m * math.comb(n + dreg - 3, dreg - 2) - first_non_trivial
    columns = R.monomials_of_degree(dreg).index(R.monomials_of_degree(dreg-2)[(rows-1)//m] * R.gens()[0] * R.gens()[0]) + 1

    return columns - rows + syzygies == planted


print(verify_monomial_semi_regularity(10, 15, planted=True))