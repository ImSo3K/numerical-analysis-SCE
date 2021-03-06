from math import exp
import numpy as np


def simpsons_one_third(f, a, b, N) -> float:
    """Calculates and prints the process of Simpson's one-third method for calculating an integral.


    Args:
        f (function): the function that is going to be integraled.
        a (int): lower boundry of the integral.
        b (int): upper boundry of the integral.
        N (int): number of subintervals of [a,b].

    Raises:
        ValueError: incase a non-even N was given.

    Returns:
        float: approximation of the integral of f(x).
    """

    if N % 2 == 1:
        raise ValueError('N must be an even integer.')

    integral = f(a)
    deltaX = (b - a) / N
    print(f'deltaX={deltaX}')
    print(f'h/3 * ({integral})', end='')
    multiplier = 2
    for i in range(1, N-1):
        fxi = multiplier*f(a+i*deltaX)
        print(f' + {multiplier}*f({a+i*deltaX:.5f}) ', end='')
        integral += fxi
        multiplier = 4 if multiplier == 2 else 2
    print(f'+ {f(b):.5f})')

    return (deltaX/3)*(integral + f(b))


def romberg(f, a, b, p) -> nd.array:
    """Calculates and prints the process of Romberg's method for calculating an integral.

    Args:
        f (function): the function that is going to be integraled.
        a (int): lower boundry of the integral.
        b (int): upper boundry of the integral.
        p (int): number of rows in the Romberg table.

        Returns:
            ndarray: 'L' shape table of approximation.
    """
    def trapezcomp(f, a, b, N) -> float:
        """composite trapezoidal function integration.
        Args:
        f (function): the function that is going to be integraled.
        a (int): lower boundry of the integral.
        b (int): upper boundry of the integral.
        N (int): number of panels in range [a,b].

        Returns:
            float: approximation of the integral of f(x).
        """

        # Initialization
        deltaX = (b - a) / N
        x = a

        # Composite rule
        In = f(a)
        for k in range(1, N):
            x = x + deltaX
            In += 2*f(x)
        return (In + f(b))*deltaX*0.5

    approx_table = np.zeros((p, p))
    for k in range(0, p):
        # Composite trapezoidal rule for 2^k panels
        approx_table[k, 0] = trapezcomp(f, a, b, 2**k)

        # Romberg recursive formula
        for j in range(0, k):
            approx_table[k, j+1] = (4**(j+1) * approx_table[k, j] - approx_table[k-1, j]) / (4**(j+1) - 1)

        print(approx_table[k, 0:k+1]) # display intermediate results

    return approx_table


def func(x): return x**2*exp(-x**2-5*x-3)*(3*x-1)


if __name__ == '__main__':
    print(f'Simpsons Solution: {simpsons_one_third(func, 0.5, 1, 20)}')
    print()
    p_rows = 2
    final_I = romberg(func, 0.5, 1, p_rows)
    solution = final_I[p_rows-1, p_rows-1]
    print(f'Romberg Solution: {solution}')
