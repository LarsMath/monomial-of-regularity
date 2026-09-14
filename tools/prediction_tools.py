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

    # Corollary 3
    if d < degree_of_regularity(n, m-1) + 2:
        return hilbert(n, m, d, start=2)
    else:
        # Corollary 4
        return macaulay_rank(n, m-1, d-2) + trivial_syzygies(n, m-1, d)
    
# ===========================================================================================================

def a_sequences(n, d):
    if d == 0: 
        yield tuple()
    else:
        for i in range(n, 0, -1):
            for a in a_sequences(i, d-1):
                yield a + (i,)

def a_sequence(n, d, i):
    assert(i <= math.comb(n+d-1,d) and i > 0)    
    remainder = math.comb(n+d-1,d) - i
    results = tuple()
    i = n
    while len(results) < d and i >= 0:
        if math.comb(i+d-len(results)-1,d-len(results)) <= remainder:
            remainder -= math.comb(i+d-len(results)-1,d-len(results))
            results = (i+1,) + results
        else:
            i -= 1
    return results

def a_index(a, n):
    return math.comb(n + len(a) - 1, len(a)) - sum(math.comb(a_k + (k+1) - 2, (k+1)) for k, a_k in enumerate(a))

# ===========================================================================================================

def nullity_predictions(n, m, solution=False, max_cols=math.inf, fast=False):
    dreg = degree_of_regularity(n, m, solution)

    while math.comb(n + dreg - 2, dreg-1) > max_cols:
        dreg -= 1

    # The lower bound is tight here
    ksyz_d = trivial_syzygies(n, m, dreg)

    predictions = {}
    # By default, go through all monomials in grevlex order
    # If fast, the variable factors of the monomial are determined in
    # sequence, from most significant to least significant
    if fast and (dreg > 2):
        tail_a = (n,)
        while True:
            a = tuple([ 1 for _ in range(dreg-2-len(tail_a))]) + tail_a
            columns = a_index((1,1) + a, n)
            i = a_index(a,n)
            rows = m * i

            # Theorem 3
            syzygies = ksyz_d - sum(trivial_syzygies(a_k - 1, m, (k+1) + 2) for k, a_k in enumerate(a))

            # Corollary 2
            if columns - rows + syzygies <= (1 if solution else 0) or (columns > max_cols): 
                if len(tail_a) == dreg-2:
                    predictions[i] = (columns, rows, syzygies)
                    break
                else: tail_a = (tail_a[0],) + tail_a
            else:
                tail_a = (tail_a[0]-1,) + tail_a[1:]
                predictions[i+1] = (columns, rows, syzygies)
    else:
        for i, a in enumerate(a_sequences(n, dreg - 2)):
            columns = a_index((1,1) + a, n)
            rows = m * (i+1)
            # Theorem 3
            syzygies = ksyz_d - sum(trivial_syzygies(a_k - 1, m, (k+1) + 2) for k, a_k in enumerate(a))

            predictions[i+1] = (columns, rows, syzygies)

            # Corollary 2
            if columns - rows + syzygies <= (1 if solution else 0) or columns > max_cols: 
                break

    return predictions

def matrix_size_at_monomial_of_regularity(n, m, solution=False):
    predictions = nullity_predictions(n, m, solution=solution)
    columns, rows, _ =  predictions[max(predictions.keys())]
    return columns, rows