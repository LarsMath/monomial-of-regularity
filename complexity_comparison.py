import math
from prediction_tools import *

def XL_complexity(n, m):
    dreg = degree_of_regularity(n, m, solution=True)    
    return 3 * math.comb(n + 1, 2) * math.comb(n + dreg - 1, dreg)**2, dreg

def FXL_complexity(n, m, q):
    cost = q**n
    best_hypers = (n, 0)
    for k in range(n):
        dreg = degree_of_regularity(n-k, m, solution=True)
        new = q**k * 3 * math.comb(n - k + 1, 2) * math.comb(n - k + dreg - 1, dreg)**2
        if new <= cost:
            best_hypers = (k, dreg)
            cost = new

    k, d = best_hypers
    return cost, k, d

def XS_complexity(n, m, q):
    cost = q**n
    best_hypers = (n, 0)
    for k in range(n-1, -1, -1):
        wiedemann = q**k * 3 * math.comb(n - k + 1, 2)
        max_cols = math.sqrt(cost / wiedemann)

        dreg = degree_of_regularity(n - k, m, solution=True)

        # For this estimation we use the heuristic that deg(mu_reg) == d_reg
        if math.comb(n - k + dreg - 2, dreg - 1) > max_cols: 
            continue

        predictions = nullity_predictions(n - k, m, solution=True, max_cols=int(max_cols))

        s = max(predictions.keys())
        c, r, syz = predictions[s]

        if c - r + syz > 1: continue

        new = wiedemann * c**2

        if new <= cost:
            best_hypers = (k, s)
            cost = new

    k, s = best_hypers
    return cost, k, s

# ===========================================================================================================================

def comparison(cases):
    results = []
    for n, m, q in cases:
        cost_xl, d_xl = XL_complexity(n, m)
        cost_fxl, k_fxl, d_fxl = FXL_complexity(n, m, q)
        cost_xs, k_xs, s_xs = XS_complexity(n, m, q)
        results += [(n, m, q, cost_xl, d_xl, cost_fxl, k_fxl, d_fxl, cost_xs, k_xs, s_xs)]

        print(f"{n} & {m} & {q} & {math.log2(cost_xl):.1f} & {d_xl} & {math.log2(cost_fxl):.1f} & {k_fxl} & {d_fxl} & {math.log2(cost_xs):.1f} & {k_xs} & {s_xs} \\\\")
    return results

q = 31
cases = [(n, n, q) for n in range(8, 31)]

results = comparison(cases)