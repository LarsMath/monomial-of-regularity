import math
from tools.prediction_tools import *
load("./tools/algebraic_tools.sage")


# In this method we verify that random systems are monomial semi-regular
# For systems with a solution we verify that they are mu-regular up to nullity 1
# We do this by finding the lowest mu with a non-trivial kernel vector
# At this point, the nullity should be at most binom(o,d)
def verify_monomial_semi_regularity(n, m, d=False, o=0, q=256):
    # Note that we dehomogenize!
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order='degrevlex')
    d_reg = degree_of_regularity(n, m, oildim=o)
    if not d:
        d = d_reg
    # Construct system
    fs, _ = construct_quadratic_system(R, n, m, o=o)

    # Construct Macaulay matrix
    mac = construct_Macaulay_matrix(R, d, fs)

    # Pivots of the actual syzygies
    ker_mac_pivots = mac.transpose().right_kernel_matrix().pivots()

    if d >= 4:
        rel = construct_relation_matrix(R, d, fs)
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
    rows = m * math.comb(n, d - 2) - first_non_trivial_pivot
    row_mon = a_sequence(n, d-2, ceil(rows/m))
    col_mon = a_shadow(row_mon,2)
    columns = a_index(col_mon,n)
    fixed_var_count = len([j for j, a_j in enumerate(row_mon) if a_j == n-d+3+j])

    return columns - rows + syzygies <= binomial(o,d)

print(verify_monomial_semi_regularity(10, 15))