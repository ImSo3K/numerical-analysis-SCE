import numpy as np


def dominant_diagonal(matrix) -> bool:
    D = np.diag(np.abs(matrix))  # Find diagonal coefficients
    S = np.sum(np.abs(matrix), axis=1) - D  # Find row sum without diagonal
    if np.all(D > S):
        print('matrix is diagonally dominant')
        return True
    else:
        print('NOT diagonally dominant')
    return False


def gauss_seidel(matrix, RHS_vec) -> None:
    if not dominant_diagonal(matrix):
        return

    ITERATION_LIMIT = 1000

    print("System of equations:")
    for i in range(matrix.shape[0]):
        row = ["{0:3g}*x{1}".format(matrix[i, j], j + 1)
               for j in range(matrix.shape[1])]
        print("[{0}] = [{1:3g}]".format(" + ".join(row), RHS_vec[i]))

    x = np.zeros_like(RHS_vec)
    for it_count in range(1, ITERATION_LIMIT):
        x_new = np.zeros_like(x)
        for i in range(matrix.shape[0]):
            s1 = np.dot(matrix[i, :i], x_new[:i])
            s2 = np.dot(matrix[i, i + 1:], x[i + 1:])
            x_new[i] = (RHS_vec[i] - s1 - s2) / matrix[i, i]
        if np.allclose(x, x_new, rtol=1e-3):  # 0.001 iter limit
            break
        x = x_new
    print("Solution: {0}".format(x))


def jacobi(matrix, RHS_vec, epsilon=1e-3, max_iterations=1000) -> None:
    if not dominant_diagonal(matrix):
        return


    print("System of equations:")
    for i in range(matrix.shape[0]):
        row = ["{0:3g}*x{1}".format(matrix[i, j], j + 1)
                for j in range(matrix.shape[1])]
        print("[{0}] = [{1:3g}]".format(" + ".join(row), RHS_vec[i]))

    x = np.zeros_like(RHS_vec)
    D = np.diag(np.diag(matrix))
    LU = matrix - D
    D_inv = np.diag(1 / np.diag(D))
    for i in range(max_iterations):
        x_new = np.dot(D_inv, RHS_vec - np.dot(LU, x))
        x = x_new
    print("Solution: {0}".format(x))


# initialize the matrix
A1 = np.array([[10., -1., 2., 0.],
               [-1., 11., -1., 3.],
               [2., -1., 10., -1.],
               [0., 3., -1., 8.]])

A2 = np.array([[3., -2., 1.],
               [1., -3., 2.],
               [-1., 2., 4.]])

# initialize the RHS vector
b1 = np.array([6.0, 25.0, -11.0, 15.0])

b2 = np.array([6.0, 25.0, -11.0])

gauss_seidel(A1, b1)
print()
gauss_seidel(A2, b2)
print()
jacobi(A1, b1)
print()
jacobi(A2, b2)
