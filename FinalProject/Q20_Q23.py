import numpy as np


def dominant_diagonal(matrix) -> bool:
    """checks if the given matrix has a dominant diagonal or not.
    Args:
        matrix (ndarray): our NxN coefficients matrix.
    Returns:
        bool: True if the matrix has a dominant diagonal, False otherwise.
    """

    D = np.diag(np.abs(matrix))  # Find diagonal coefficients
    S = np.sum(np.abs(matrix), axis=1) - D  # Find row sum without diagonal
    if np.all(D > S):
        print('matrix is diagonally dominant')
        return True
    else:
        print('NOT diagonally dominant')
    return False


def gauss_seidel(matrix, RHS_vec) -> None:
    """solves a system of linear equations with Gauss Seidel method.
    Args:
        matrix (ndarray): our coefficients NxN matrix.
        RHS_vec (ndarray): solution vector {b}.
    """
    if not dominant_diagonal(matrix):
        return

    ITERATION_LIMIT = 1000

    print('System of equations:')
    for i in range(matrix.shape[0]):
        row = ['{0:3g}*x{1}'.format(matrix[i, j], j + 1)
               for j in range(matrix.shape[1])]
        print('[{0}] = [{1:3g}]'.format(' + '.join(row), RHS_vec[i]))
    iter_no = 0
    x = np.zeros_like(RHS_vec)
    for it_count in range(1, ITERATION_LIMIT):
        x_new = np.zeros_like(x)
        print('iterative equetions:')
        for i in range(matrix.shape[0]):
            s1 = np.dot(matrix[i, :i], x_new[:i])
            s2 = np.dot(matrix[i, i + 1:], x[i + 1:])
            x_new[i] = (RHS_vec[i] - s1 - s2) / matrix[i, i]
            print(f'[x{i+1}r+1 = ({RHS_vec[i]} - {s1}*x2 - {s2}*x3)/{matrix[i,i]}]')
        if np.allclose(x, x_new, rtol=1e-5):  # 0.001 iter limit
            break
        x = x_new
        print(
            f'[x1r+1 = {x[0]}, x2r+1 = {x[1]}, x3r+1 = {x[2]}] iteration No. : {it_count}')
        iter_no = it_count
    print(
        f'Solution found : [x1r+1 = {x[0]}, x2r+1 = {x[1]}, x3r+1 = {x[2]}] iteration No. : {iter_no}')
    return x


def jacobi(matrix, RHS_vec, epsilon=1e-5, max_iterations=1000) -> None:
    """solves a system of linear equations with Jacobi method.
    Args:
        matrix (ndarray): our coefficients NxN matrix.
        RHS_vec (ndarray): solution vector {b}.
    """
    if not dominant_diagonal(matrix):
        return

    print('System of equations:')
    for i in range(matrix.shape[0]):
        row = ['{0:3g}*x{1}'.format(matrix[i, j], j + 1)
               for j in range(matrix.shape[1])]
        print('[{0}] = [{1:3g}]'.format(' + '.join(row), RHS_vec[i]))
    iter_no = 0
    x = np.zeros_like(RHS_vec)
    D = np.diag(np.diag(matrix))
    LU = matrix - D
    D_inv = np.diag(1 / np.diag(D))
    for i in range(max_iterations):
        x_new = np.dot(D_inv, RHS_vec - np.dot(LU, x))
        if np.allclose(x, x_new, rtol=1e-5):
            break
        x = x_new
        print(
            f'[x1r+1 = {x[0]}, x2r+1 = {x[1]}, x3r+1 = {x[2]}] iteration No. : {i}')
        iter_no = i
    print(
        f'Solution found : [x1r+1 = {x[0]}, x2r+1 = {x[1]}, x3r+1 = {x[2]}] iteration No. : {iter_no}')

# initialize the matrix


# Question No. 20:
A3 = np.array([[10, 8, 1],
               [4, 10, -5],
               [5, 1, 10]])
# Question No. 23:
A2 = np.array([[0.04, 0.01, 0.01],
               [0.2, 0.5, -0.2],
               [1, 2, 4]])

# initialize the RHS vector
# Question No. 20:
b3 = np.array([-7, 2, 1.5])
# Question No. 23:
b2 = np.array([0.06, 0.3, 11])

if __name__ == '__main__':
    print('#Question No. 23 using gauss_seidel:')
    gauss_seidel(A2, b2)  # 23
    print()
    print('#Question No. 23 using jacobi:')
    jacobi(A2, b2)
    print()
    print('#Question No. 20 using gauss_seidel:')
    gauss_seidel(A3, b3)  # 20
    print()
    print('#Question No. 20 using jacobi:')
    jacobi(A3, b3)
