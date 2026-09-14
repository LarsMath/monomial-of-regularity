import math
from functools import cache

# ===================================== Combinatorial homological algebra =================================================

@cache
def hilbert(n, m, d, start=0):
    return sum((-1)**i * (math.comb(m+i-1, i) if m>0 else 0 if i>0 else 1) * math.comb(n, d - 2*i) for i in range(start, d//2 + 1))

@cache
def degree_of_regularity(n, m, oildim = 0):
    for d in range(2,n+1):
        if hilbert(n, m, d) <= math.comb(oildim,d):
            return d
    return math.inf

@cache
def macaulay_nullity(n, m, d):
    if d < 2: return 0
    if d > 2 and macaulay_nullity(n, m, d-1) == 0: return 0
    
    return max(0, hilbert(n, m, d))

@cache
def macaulay_rank(n, m, d):
    if m <= 0: return 0
    return math.comb(n, d) - macaulay_nullity(n, m, d)
    
# Note that this is in general a lower bound!
@cache
def trivial_syzygies(n, m, d):
    if d < 4 or m <= 0: return 0

    # Corollary 6
    return max(macaulay_rank(n, m, d-2),0) + max(trivial_syzygies(n, m-1, d),0)
    
# ===========================================================================================================

def a_sequences(n, d):
    if d == 0: yield tuple()
    for i in range(n, 0, -1):
        for a in a_sequences(i-1, d-1):
            yield a + (i,)

def a_sequence(n, d, i):
    assert(i <= math.comb(n,d) and i > 0)    
    remainder = math.comb(n,d) - i
    results = tuple()
    i = n
    while len(results) < d and i >= 0:
        if math.comb(i,d-len(results)) <= remainder:
            remainder -= math.comb(i,d-len(results))
            results = (i+1,) + results
        i -= 1
    return results

def a_index(a, n):
    return math.comb(n, len(a)) - sum(math.comb(a_k-1, (k+1)) for k, a_k in enumerate(a))

def a_shadow(a, c):
    shadow_a = a
    i = 0
    while len(shadow_a) < len(a) + c:
        if i == len(shadow_a): return tuple(range(1,len(a)+c+1))
        if shadow_a[i] != i+1: shadow_a = shadow_a[:i] + (i+1,) + shadow_a[i:]
        i+=1
    return shadow_a

# ===========================================================================================================

def nullity_predictions(n, m, d, solution=False, max_cols = math.inf, fast = False):
    # The lower bound is tight here
    ksyz_d = trivial_syzygies(n, m, d)

    predictions = {}
    
    # By default, go through all monomials in grevlex order
    # If fast, the variable factors of the monomial are determined in
    # sequence, from most significant to least significant
    if fast:
        tail_a = (n,)
        while True:
            a = tuple(range(1,d-1-len(tail_a))) + tail_a

            columns = a_index(a_shadow(a,2), n)
            i = a_index(a,n)
            rows = m * i
            fixed_var_count = len([j for j, a_j in enumerate(a) if a_j == n-d+3+j])

            # Corollary 5
            syzygies = ksyz_d - sum(trivial_syzygies(a_j - 1, m, (j+1) + 2) for j, a_j in enumerate(a))
            if columns - rows + syzygies <= (1 if solution else 0) or columns > max_cols: 
                if len(tail_a) == d-2:
                    predictions[i] = (columns, rows, syzygies, fixed_var_count)
                    break
                else: tail_a = (tail_a[0]-1,) + tail_a
            else:
                tail_a = (tail_a[0]-1,) + tail_a[1:]
                predictions[i] = (columns, rows, syzygies, fixed_var_count)
    else:
        for i, a in enumerate(a_sequences(n, d-2)):
            columns = a_index(a_shadow(a,2), n)
            rows = m * (i+1)
       
            fixed_var_count = len([j for j, a_j in enumerate(a) if a_j == n-d+3+j])

            # Corollary 5
            syzygies = ksyz_d - sum(trivial_syzygies(a_j - 1, m, (j+1) + 2) for j, a_j in enumerate(a))

            predictions[i+1] = (columns, rows, syzygies, fixed_var_count)
            if columns - rows + syzygies <= (1 if solution else 0) or columns > max_cols: break
    return predictions


def matrix_size_at_monomial_of_regularity(n, m, solution=False):
    predictions = nullity_predictions(n, m, solution=solution)
    columns, rows, _ =  predictions[max(predictions.keys())]
    return columns, rows