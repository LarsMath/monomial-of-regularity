from tools.prediction_tools import *
load("./tools/algebraic_tools.sage")

# This method computes the nullities of random quadratic systems for Mac^mu for each choice of mu in the Macaulay matrix of degree d
# Returns a dictionary with |S|: (columns, rows, syzygies)
def compute_nullities(n, m, d, q=256):
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order='degrevlex')

    fs, _ = construct_quadratic_system(R, n, m)

    mac = construct_Macaulay_matrix(R, d, fs)

    syz_pivots = mac.transpose().right_kernel_matrix().pivots()

    return {i+1: sum(1 for s in syz_pivots if m * binomial(n, d - 2) - s <= m * (i+1)) for i in range(binomial(n,d-2))}

n = 10
m = 3
d = degree_of_regularity(n,m)

prediction = nullity_predictions(n, m, 5, fast = False)
actual = compute_nullities(n, m, 5)

print("|S|\tpred\tactual\tupper bound")
for s, (columns, rows, pred_syzygies, _) in prediction.items():
    syzygies = actual[s]
    print(f"{s}\t{columns - rows + pred_syzygies}\t{columns - rows + syzygies}\t{columns - rows + pred_syzygies >= columns - rows + syzygies}")