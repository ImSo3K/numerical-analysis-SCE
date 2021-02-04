import math
from numpy import arange


def newton_raphson(f, x0, epsilon=1e-5) -> None:
    """find root of a real function with Newton-Raphson method.
    Args:
        f (function): the function whose root we are trying to find.
        x0 (float): the initial guess.
        epsilon (float, optional): tolerance value. Defaults to 1e-4.
    """
    def func_der(x):
        return (-6*math.exp(math.pow(-x, 2)-5*x-3)*math.pow(x, 4)-13*math.exp(math.pow(-x, 2)-5*x-3)*math.pow(x, 3)+14*math.exp(math.pow(-x, 2)-5*x-3)*math.pow(x, 2)-2*math.exp(math.pow(-x, 2)-5*x-3)*x)

    x1 = 0
    solution = False
    for i in range(100):
        y = f(x0)
        y_der = func_der(x0)

        x1 = x0 - y/y_der
        print(f'iteration No. {i+1}, x1 = {x0} - {y/y_der}')
        if abs(x1 - x0) <= epsilon:
            solution = True
            print(f'solution found : {x1}, iteration number : {i+1}')
            return

        if y_der == 0:
            print('Zero derivative. No solution found.')
            return

        x0 = x1
    if not solution:
        print('Did not converge')


def secant(f, x0, x1, epsilon=1e-5) -> None:
    """find root of a real function with Newton-Raphson method.
    Args:
        f (function): the function whose root we are trying to find.
        x0 (float): the first initial guess.
        x1 (float): the second initial guess.
        epsilon (float, optional): tolerance value. Defaults to 1e-4.
    """
    for i in range(0, 100):
        x = x1 - ((f(x1)*(x1-x0))/(f(x1) - f(x0)))
        print(f'iteration No. {i+1}, x = {x1} - {((f(x1)*(x1-x0))/(f(x1) - f(x0)))}')
        if abs(x-x1) < epsilon:
            print(f'solution found : {x}, iteration number : {i+1}')
            return
        x0 = x1
        x1 = x


def func(x): return (math.pow(x, 2)*math.exp(math.pow(-x, 2)-5*x-3))*(3*x-1)


def activator(method, f) -> None:
    """activates the proper method to find the roots of a function.
    Args:
        method (function): which method we want to use.
        f (function): which function is going to be used with the method.
    """
    item = range_tuple[0]
    stp = 0.01
    for i in arange(range_tuple[0], range_tuple[1], stp):
        if (f(item) > 0 and f(item+stp) < 0) or (f(item) < 0 and f(item+stp) > 0):
            print(f'roots in range [{item}, {item + stp}] : ')
            if method.__name__ == newton_raphson.__name__:
                method(f, item + stp, 1e-5)
            else:
                method(f, item, item + stp, 1e-5)
        item += stp
    if(func(0) == 0):
        print(f'solution found : {func(0)}, iteration number : {0}')


if __name__ == "__main__":
    range_tuple = (0, 1.5)
    print('Newton-Rhapson roots: ')
    activator(newton_raphson, func)
    print()
    print('Secant roots: ')
    activator(secant, func)
