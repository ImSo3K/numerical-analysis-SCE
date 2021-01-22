import numpy as np
from Q20_Q23 import gauss_seidel


def cubic_spline(data, value) -> None:
    """cubic splining to interpolate a value according to given data points.

    Args:
        data (tuple): tuple of tuples which formated as '(x_coor, y_coor)'.
        value (float): the value to calculate it's interpolation.
    Returns:
        list: a list containing the results of each function corresponding to each pair of points
    """
    if len(data) < 2:
        return 'failed'

    def x(i) -> float:
        return data[i][0]

    def y(i) -> float:
        return data[i][1]

    def h(i) -> float:
        return x(i + 1) - x(i)

    n = len(data)
    matrix_A = np.identity(n)
    sol_vector_b = [(y(i + 1) - y(i)) / (h(i)) - (y(i) - y(i - 1)) / (h(i - 1))
                    for i in range(1, n - 1)]

    sol_vector_b.insert(0, 0)
    sol_vector_b.append(0)
    for i in range(1, n - 1):
        matrix_A[i][i - 1] = h(i - 1) / 6
        matrix_A[i][i] = ((h(i - 1) + h(i)) / 3)
        matrix_A[i][i + 1] = h(i) / 6

    print(
        f'The matrix:\n{matrix_A}\n\nThe solution vector of the cubic spline:\n{sol_vector_b}\n')
    print('calculating the matrix results with gauss seidel')

    Mi = gauss_seidel(matrix_A, sol_vector_b)
    for i, j in enumerate(Mi):
        print(f'M({i}) = {j}')
    print()
    results = [0.0 for _ in range(n - 1)]
    for i in range(n - 1):
        print(
            f'S{i}(x)= ({y(i + 1)}(x-{x(i)})/{h(i)})-({y(i)}(x-{x(i + 1)})/{h(i)}) + ({Mi[i + 1]}/6)[((x-{x(i)}) ^ 3)/{h(i)}) - {h(i)}(x-{x(i)}))] -({Mi[i]}/6)[((x-{x(i + 1)}) ^ 3)/{h(i)}) - {h(i)}(x-{x(i + 1)})]\n')

        results[i] = (y(i + 1) * (value - x(i)) / h(i)) - (y(i) * (value - x(i + 1)) / h(i)) + \
            (Mi[i + 1] / 6) * ((((value - x(i)) ** 3) / h(i)) - h(i) * (value - x(i))) - \
            (Mi[i] / 6) * ((((value - x(i + 1)) ** 3) / h(i)) -
                           h(i) * (value - x(i + 1)))
        print(f'S{i}({value}) = {results[i]}\n')


if __name__ == '__main__':
    cubic_spline(((0.1, -0.29004996), (0.2, -0.56079734),
                  (0.3, -0.81401972)), 0.25)
