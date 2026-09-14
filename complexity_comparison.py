import math
from tools.prediction_tools import *

def XL_complexity(n, m):
    dreg = degree_of_regularity(n, m, solution=False)    
    return 3 * math.comb(n + 1, 2) * math.comb(n + dreg - 1, dreg)**2, dreg

def FXL_complexity(n, m, q):
    cost = q**n
    best_hypers = (n, 0)
    for k in range(max(0,n-m),n):
        dreg = degree_of_regularity(n-k, m, solution=False)
        new = q**k * 3 * math.comb(n - k + 1, 2) * math.comb(n - k + dreg - 1, dreg)**2
        if new <= cost:
            best_hypers = (k, dreg)
            cost = new

    k, d = best_hypers
    return cost, k, d

def XS_complexity(n, m, fast=False):
    wiedemann = 3 * math.comb(n + 1, 2)

    predictions = nullity_predictions(n, m, solution=False, fast=fast)

    s = max(predictions.keys())
    c, r, syz = predictions[s]

    if c - r + syz > 1: return math.inf

    cost = wiedemann * c**2
    return cost, s

def FXS_complexity(n, m, q, fast=False):
    cost = q**n
    best_hypers = (n, 0)
    for k in range(n-1, max(0,n-m)-1, -1):
        wiedemann = q**k * 3 * math.comb(n - k + 1, 2)
        max_cols = math.sqrt(cost / wiedemann)

        dreg = degree_of_regularity(n - k, m, solution=True)

        # For this estimation we use the heuristic that deg(mu_reg) == d_reg
        if math.comb(n - k + dreg - 2, dreg - 1) > max_cols: 
            continue

        predictions = nullity_predictions(n - k, m, solution=True, max_cols=int(max_cols), fast=fast)

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

def comparison(cases, fast=False):
    results = []
    for n, m, q in cases:
        cost_xl, d_xl = XL_complexity(n, m)
        cost_xs, s_xs = XS_complexity(n, m, fast=fast)
        results += [(n, m, q, cost_xl, d_xl, cost_xs, s_xs)]

        print(f"{n} & {m} & {q} & {d_xl} & {math.log2(cost_xl):.1f} & {math.log2(cost_xs):.1f} \\\\")
    return results

def hybrid_comparison(cases, fast=False):
    results = []
    for n, m, q in cases:
        q_factor = math.log2(q) * q
        cost_fxl, k_fxl, d_fxl = FXL_complexity(n, m, q)
        cost_fxs, k_fxs, s_fxs = FXS_complexity(n, m, q, fast=fast)
        results += [(n, m, q, cost_fxl, k_fxl, d_fxl, cost_fxs, k_fxs, s_fxs)]

        print(f"{n} & {m} & {q} & {math.log2(cost_fxl * q_factor):.1f} & {k_fxl} & {d_fxl} & {math.log2(cost_fxs* q_factor):.1f} & {k_fxs} & {s_fxs} \\\\")
    return results

q = 31
cases = [(n, n, q) for n in range(5, 36)]
results = hybrid_comparison(cases, fast=True)