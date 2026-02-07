from tools.prediction_tools import *
load("./tools/algebraic_tools.sage")

# This takes a quadratic monomial semi-regular system fs (with a solution) and a ring with the degrevlex order
# It then computes a solution by predicting the monomial of regularity and solving the associated Macaulay matrix
# It will fail is the nullity is not 1
def sparse_solver(R, fs):

    n, m = len(R.gens()), len(fs)

    d_reg = degree_of_regularity(n, m, solution=True)
    columns, rows = matrix_size_at_monomial_of_regularity(n, m, solution=True)

    # Construct Macaulay matrix
    mac = construct_Macaulay_matrix(R, d_reg, fs)

    # One could optimize here to not construct the entire Macaulay matrix in the first place
    mac = mac.matrix_from_rows_and_columns([i for i in range(mac.nrows() - rows, mac.nrows())], [i for i in range(mac.ncols() - columns, mac.ncols())])

    found = mac.right_kernel()

    assert found.dimension() == 1, "The right nullity is higher than 1! This system is not monomial semi-regular"

    found_solution = found.basis()[0]
    found_solution = [x / found_solution[-1] for x in found_solution[-n:]]
    return found_solution

n = 7
m = 7
R = PolynomialRing(GF(251), n, 'x', order='degrevlex')

fs, sol = construct_quadratic_system(R, n, m, True)
found = sparse_solver(R, fs)

print(sol, found, sol==found) 