import math
from functools import cache

# ===================================== Combinatorial homological algebra =================================================

# Everything is homogeneous!
@cache
def hilbert(n, m, d, start=0):
    return sum((-1)**i * math.comb(m, i) * math.comb(n + d - 2*i - 1, d - 2*i) for i in range(start, d//2 + 1))

@cache
def degree_of_regularity(n, m, solution=False):
    return next(d for d in range(2, n+2) if hilbert(n, m, d) <= (1 if solution else 0)) if n <= m else math.inf

@cache
def macaulay_nullity(n, m, d):
    if d < 2: return 0
    if d > 2 and macaulay_nullity(n, m, d-1) == 0: return 0
    
    return max(0, hilbert(n, m, d))

@cache
def macaulay_rank(n, m, d):
    if m <= 0: return 0
    return math.comb(n + d - 1, d) - macaulay_nullity(n, m, d)
    
# Note that this is in general a lower bound!
@cache
def trivial_syzygies(n, m, d):
    if d < 4 or m <= 0 or n <= 0: return 0

    # Corollary 4.18
    if d < degree_of_regularity(n, m-1) + 2:
        return hilbert(n, m, d, start=2)
    else:
        # Corollary 4.16
        return macaulay_rank(n, m-1, d-2) + trivial_syzygies(n, m-1, d)
    
# ===========================================================================================================

@cache
def a_sequences(n, d):
    if d == 0: return [tuple()]
    return [a + (i, ) for i in range(n, 0, -1) for a in a_sequences(i, d-1)]

@cache
def a_index(a, n):
    return math.comb(n + len(a) - 1, len(a)) - sum(math.comb(a_k + (k+1) - 2, (k+1)) for k, a_k in enumerate(a))

# ===========================================================================================================

def nullity_predictions(n, m, solution=False, max_cols=math.inf):
    dreg = degree_of_regularity(n, m, solution)

    while math.comb(n + dreg - 2, dreg-1) > max_cols:
        dreg -= 1

    # The lower bound is tight here
    ksyz_d = trivial_syzygies(n, m, dreg)

    predictions = {}
    for i, a in enumerate(a_sequences(n, dreg - 2)):

        columns = a_index((1,1) + a, n)
        rows = m * (i+1)
        # Theorem 4.13
        syzygies = ksyz_d - sum(trivial_syzygies(a_k - 1, m, (k+1) + 2) for k, a_k in enumerate(a))

        predictions[i+1] = (columns, rows, syzygies)

        # Corollary 4.14
        if columns - rows + syzygies <= (1 if solution else 0) or columns > max_cols: break

    return predictions