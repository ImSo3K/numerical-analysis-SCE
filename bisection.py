from math import log, ceil
import numpy as np
import itertools


def bisection(start_point, end_point, f, epsilon):
    """finding the root of the function with bisection method.
    Args:
        start_point (float): intuitive starting point.
        end_point (float): intuitive end point.
        f (function): the function which we will try to find it's root.
        epsilon (float): tolerance
    """
    if f(start_point) * f(end_point) > 0.0:
        print('start point & end value are invalid, try again with different guess values')
        return
    m = (start_point + end_point) / 2.0
    # number of iterations upper bound
    k = ceil(-(log(epsilon/(end_point-start_point))/log(2)))
    for i in range(k):

        if abs(start_point - end_point) > epsilon:
            if f(start_point) * f(m) < 0:
                end_point = m
            else:
                start_point = m
            m = (start_point + end_point) / 2.0
            print(
                f"iteration number: {i+1}: {start_point} + {end_point} / 2={m}")

    # rounding works until 1e-4 because from there the number repr changes to "1e-x"
    if(f==f_derivative and function(m)==0):
        roots.append((round(m,len(str(epsilon)[2:]))))
    elif abs(f(round(m, int(str(epsilon).split('.')[1])))) == 0 or abs(start_point - end_point) < epsilon and f!=f_derivative:
        roots.append((round(m,len(str(epsilon)[2:]))))
    else:
        print('no root for this function in given')


def range_bisection(min_range, max_range, section):
    """send each valid section to the bisection algorithm.
    Args:
        min_range (int): starting point of the total range
        max_range (int): ending point of the total range
        section (float): the gap between each pair of sub-sections.
    """
    number_of_sections = (max_range - min_range) / section
    section_list = np.linspace(min_range, max_range, int(number_of_sections)+1)
    # section_list = section_list[section_list != 0]
    a, b = itertools.tee(section_list)
    next(b, None)
    pair = zip(a, b)
    for ranges in pair:
        if function(ranges[0]) * function(ranges[1]) < 0.:
            bisection(ranges[0], ranges[1], function, epsilon)
        if f_derivative(ranges[0]) * f_derivative(ranges[1]) < 0.:
            bisection(ranges[0], ranges[1], f_derivative, epsilon)
    if (function(0)) == 0:
        roots.append(0)


def function(x): return x**4 + x**3 - 3*x**2
def f_derivative(x): return 4*x**3 + 3*x**2 - 6*x


if __name__ == "__main__":
    epsilon = 1e-4
    roots = []
    range_bisection(-3.0, 2.0, 0.1)
    print(f'roots : {roots}')
    roots.clear()