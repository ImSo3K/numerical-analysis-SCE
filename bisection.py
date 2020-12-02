from math import log, ceil
def bisection(start_point, end_point,f, epsilon) -> float:
    """finding the root of the function with bisection method.
    Args:
        start_point (float): intuitive starting point.
        end_point (float): intuitive end point.
        f (function): the function which we will try to find it's root.
        epsilon (float): tolerance

    Returns:
        float: the function's root.
    """
    m = (start_point + end_point) / 2.0
    k = ceil(-(log(epsilon/(end_point-start_point))/log(2))) # number of iterations upper bound
    for i in range(k):
        if abs(start_point - end_point) > epsilon:
            if f(m) == 0:
                print(f"took {i} iterations")
                return m
            elif f(start_point) * f(m) < 0:
                end_point = m
            else:
                start_point = m
            m = (start_point + end_point) / 2.0
            # print(f"{start_point} + {end_point} / 2={m}")
        # print(f"took {i} iterations")
    return m

print(bisection(-1.5,1.5,lambda x: x*x +3*x + 2, 1e-4))