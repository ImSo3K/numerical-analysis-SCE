from math import cos


def newton_raphson(f, x0, epsilon=1e-4):
    """find root of a real function with Newton-Raphson method.

    Args:
        f (function): the function whose root we are trying to find.
        x0 (float): the initial guess.
        epsilon (float, optional): tolerance value. Defaults to 1e-4.
    """
    def func_der(x):
        return 12*x**2 - 48

    x1 = 0
    solution = False
    for i in range(100):
        y = f(x0)
        y_der = func_der(x0)

        x1 = x0 - y/y_der
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


def secant(f, x0, x1, epsilon=1e-4):
    """find root of a real function with Newton-Raphson method.

    Args:
        f (function): the function whose root we are trying to find.
        x0 (float): the first initial guess.
        x1 (float): the second initial guess.
        epsilon (float, optional): tolerance value. Defaults to 1e-4.
    """
    for i in range(0, 100):
        x = x1 - ((f(x1)*(x1-x0))/(f(x1) - f(x0)))
        if abs(x-x1) < epsilon:
            print(f'solution found : {x}, iteration number : {i+1}')
            return
        x0 = x1
        x1 = x


def func(x): return 4*x**3 - 48*x + 5


def activator(method, f):
    """activates the proper method to find the roots of a function.

    Args:
        method (function): which method we want to use.
        f (function): which function is going to be used with the method.
    """
    item = range_tuple[0]
    if range_tuple[1] - range_tuple[0] >= 1:
        for i in range(range_tuple[1]):
            if (f(item) > 0 and f(item+1) < 0) or (f(item) < 0 and f(item+1) > 0):
                print(f'roots in range [{item}, {item + 1}] : ')
                if method.__name__ == newton_raphson.__name__:
                    method(f, item + 0.5, 1e-4)
                else:
                    method(f, item, item + 1, 1e-3)
            item += 1


if __name__ == "__main__":
    range_tuple = (0, 10)
    print('Newton-Rhapson roots: ')
    activator(newton_raphson, func)
    print()
    print('Secant roots: ')
    activator(secant, lambda x: x**3 - cos(x))
