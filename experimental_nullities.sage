from tools.prediction_tools import *
load("./tools/algebraic_tools.sage")

# This method computes the nullities of random quadratic systems for Mac^mu for each choice of mu in the Macaulay matrix of degree d
# Returns a dictionary with |S|: (columns, rows, syzygies)
def compute_nullities(n, m, d, q=251):
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order='degrevlex')

    fs, _ = construct_quadratic_system(R, n, m)

    mac = construct_Macaulay_matrix(R, d, fs)

    syz_pivots = mac.transpose().right_kernel_matrix().pivots()

    return {i+1: (binomial(n + d - 1, d) - R.monomials_of_degree(d)[::-1].index(mult * R.gens()[0]**2), m * (i+1), sum(1 for s in syz_pivots if m * binomial(n + d - 3, d - 2) - s <= m * (i+1))) for i, mult in enumerate(R.monomials_of_degree(d-2))}

n = 7
m = 8

prediction = nullity_predictions(n, m)
actual = compute_nullities(n, m, 5)

print("|S|\tpred\tactual\tupper bound")
for s, (pred_columns, pred_rows, pred_syzygies) in prediction.items():
    columns, rows, syzygies = actual[s]
    print(f"{s}\t{pred_columns - pred_rows + pred_syzygies}\t{columns - rows + syzygies}\t{pred_columns - pred_rows + pred_syzygies >= columns - rows + syzygies}")