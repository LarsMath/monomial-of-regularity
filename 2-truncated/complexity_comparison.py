import math
from tools.prediction_tools import *

def truncated_XL_complexity(n, m, o, q):
    cost = 0
    best_k = 0
    field_costs = 2* math.log(q,2)**2 + math.log(q,2)
    for k in range(o+1):
        dreg = degree_of_regularity(n-k, m, oildim = o-k)
        if dreg <= o-k:
            new = 3 * math.comb(n - k+1, 2) * math.comb(n - k, dreg)**2 * field_costs
            if new <= cost or cost == 0:
                best_k = k
                cost = new
    
    return cost, best_k

def truncated_XS_complexity(n, m, o, q, fast=False):
    wiedemann = 3 * math.comb(n+1, 2)
    field_costs = 2* math.log(q,2)**2 + math.log(q,2)
    predictions = nullity_predictions(n, m, o, solution=True, fast=fast)

    s = max(predictions.keys())
    c, r, syz, fixed = predictions[s]

    if c - r + syz > 1: return 0

    cost = wiedemann * c**2 * field_costs

    return cost, fixed, s

# ===========================================================================================================================

def comparison(cases, fast=False):
    results = []
    for n, m, o, q in cases:
        cost_txl, k_txl = truncated_XL_complexity(n, m, o,q )
        if cost_txl != 0:
            cost_txs, k_txs, s_txs = truncated_XS_complexity(n-k_txl, m, o-k_txl, q, fast=fast)
            k_txs += k_txl
            results += [(n, m, q, cost_txl, k_txl, cost_txs, k_txs, s_txs)]
            print(f"{n} & {m} & {q} & {math.log2(cost_txl):.1f} & {k_txl} & {math.log2(cost_txs):.1f} & {k_txs} & {s_txs} \\\\")
        else:
            print(f"{n} & {m} & {q} & - & - & - & - & - \\\\")
    return results


UOV = [ (112,44,44,256), (160,64,64,16), (184,72,72,256), (244,96,96,256) ]
print("UOV, Reconcilliation:")
results = comparison(UOV, fast=True)

intersection_UOV = [ (n, 3*m-2, 3*o - n, q) for n,m,o,q in UOV ]
print("UOV, Intersection:")
results = comparison(intersection_UOV, fast=True)

MAYO = [(81, 64, 17, 16)]
print("MAYO:")
results = comparison(MAYO, fast=True)