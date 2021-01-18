from math import sin, exp
import numpy as np


def simpsons_one_third(f, a, b, N):
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
        raise ValueError("N must be an even integer.")

    integral = f(a)
    deltaX = (b - a) / N
    print(f"deltaX={deltaX}")
    print(f'h/3 * ({integral})', end="")
    multiplier = 2
    for i in range(1, N-1):
        fxi = multiplier*f(a+i*deltaX)
        print(f" + {multiplier}*f({a+i*deltaX:.5f}) ", end="")
        integral += fxi
        multiplier = 4 if multiplier == 2 else 2
    print(f"+ {f(b):.5f})")

    return (deltaX/3)*(integral + f(b))


def romberg(f, a, b, p):
    """Calculates and prints the process of Romberg's method for calculating an integral.

    Args:
        f (function): the function that is going to be integraled.
        a (int): lower boundry of the integral.
        b (int): upper boundry of the integral.
        p (int): number of rows in the Romberg table.
    """
    def trapezcomp(f, a, b, N):
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

    I = np.zeros((p, p))
    for k in range(0, p):
        # Composite trapezoidal rule for 2^k panels
        I[k, 0] = trapezcomp(f, a, b, 2**k)

        # Romberg recursive formula
        for j in range(0, k):
            I[k, j+1] = (4**(j+1) * I[k, j] - I[k-1, j]) / (4**(j+1) - 1)

        print(I[k, 0:k+1])   # display intermediate results

    return I


def func(x): return (sin(2*x**3 + 5*x**2 - 6))/(2*exp(-2*x))


if __name__ == '__main__':
    p_rows = 5
    final_I = romberg(func, 0, 1, p_rows)
    solution = final_I[p_rows-1, p_rows-1]
    print(solution)
    print()
    print(simpsons_one_third(func, 0, 1, 2500))
