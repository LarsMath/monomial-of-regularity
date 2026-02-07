import math
from tools.prediction_tools import *
load("./tools/algebraic_tools.sage")


# In this method we verify that random systems are monomial semi-regular
# For systems with a solution we verify that they are mu-regular up to nullity 1
# We do this by finding the lowest mu with a non-trivial kernel vector
# At this point, the nullity should be exactly 0 or 1
def verify_monomial_semi_regularity(n, m, planted, q=251):
    # Note that we dehomogenize!
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order='degrevlex')
    d_reg = degree_of_regularity(n, m, solution=planted)

    # Construct system
    fs, _ = construct_quadratic_system(R, n, m, planted)

    # Construct Macaulay matrix
    mac = construct_Macaulay_matrix(R, d_reg, fs)

    # Pivots of the actual syzygies
    ker_mac_pivots = mac.transpose().right_kernel_matrix().pivots()

    if d_reg >= 4:
        rel = construct_relation_matrix(R, d_reg, fs)
        # Pivots of the trivial syzygies
        rel_pivots = rel.pivots()
    else:
        # No trivial syzygies
        rel_pivots = []

    if all(p in rel_pivots for p in ker_mac_pivots):
        # All syzygies are trivial
        return True
    
    first_non_trivial_pivot = next(p for p in ker_mac_pivots[::-1] if p not in rel_pivots)

    syzygies = sum(1 for p in ker_mac_pivots if p >= first_non_trivial_pivot)
    rows = m * math.comb(n + d_reg - 3, d_reg - 2) - first_non_trivial_pivot
    columns = R.monomials_of_degree(d_reg).index(R.monomials_of_degree(d_reg-2)[(rows-1)//m] * R.gens()[0] * R.gens()[0]) + 1

    return columns - rows + syzygies == planted


print(verify_monomial_semi_regularity(10, 15, planted=True))