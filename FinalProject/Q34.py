import numpy as np


def lagrange(x, y, x_for_approx):
    """interpolation for the given value 'x_for_approx'.

    Args:
        x (ndarray): numpy array which contains the x column.
        y (ndarray): numpy array which contains the f(x) (y values) column.
        x_for_approx (int): the value we are going to approximate with the interpolation.

    Returns:
        float: the approximated result.
    """
    approx_result = 0
    num_of_iteration = 0
    for i in range(0, x.size):
        poly = y[i]
        for j in range(0, x.size):
            if i != j:
                print(f'poly = {poly}*{(x_for_approx - x[j])}/{(x[i] - x[j])}')
                poly = poly * (x_for_approx - x[j]) / (x[i] - x[j])
            num_of_iteration += 1
        print(f'f`(x) = {approx_result}+{poly}, iteration No.: {num_of_iteration}')
        approx_result = approx_result + poly
        num_of_iteration += 1
    return approx_result


def neville(x, y, x_for_approx):
    """interpolation for the given value 'x_for_approx'.

    Args:
        x (ndarray): numpy array which contains the x column.
        y (ndarray): numpy array which contains the f(x) (y values) column.
        x_for_approx (int): the value we are going to approximate with the interpolation.

    Returns:
        float, ndarray: the approximated result, matrix of coefficients.
    """

    coefficients_matrix = np.zeros((x.size, x.size - 1))
    coefficients_matrix = np.concatenate((y[:, None], coefficients_matrix), axis=1)

    for i in range(1, x.size):
        for j in range(1, i + 1):
            print(
                f'coefficients_matrix[{i},{j}] = (({x_for_approx}-{x[i-j]})*{coefficients_matrix[i, j - 1]} - {(x_for_approx - x[i])} * {coefficients_matrix[i - 1, j - 1]}) / {(x[i] - x[i - j])})')
            coefficients_matrix[i, j] = ((x_for_approx - x[i - j]) * coefficients_matrix[i, j - 1] -
                       (x_for_approx - x[i]) * coefficients_matrix[i - 1, j - 1]) / (x[i] - x[i - j])

    approx_result = coefficients_matrix[x.size - 1, x.size - 1]
    return approx_result, coefficients_matrix


x_values = np.array([1.2, 1.3, 1.4, 1.5, 1.6])
y_values = np.array([3.5095, 3.6984, 3.9043, 4.1293, 4.3756])
x_for_approx = 1.37

print('Approximation in Neville`s method:')
print(f'Our x values: {x_values}')
print(f'Our x values: {y_values}')
print(f'value for approximation: {x_for_approx}')
approx_result, coefficients_matrix = neville(x_values, y_values, x_for_approx)

print(f'Coefficients matrix:\n{coefficients_matrix}')
print(f"the Neville`s interpolated value at {x_for_approx} is: {approx_result}")

print()
print('Lagrange Interpolation:')
approx_result = lagrange(x_values, y_values, x_for_approx)
print(f"the Lagrange`s interpolated value at {x_for_approx} is: {approx_result}")
