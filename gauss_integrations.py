from gauss_schemes import GaussIntegral

def function_1d(x : float) -> float:
    return 5*(x**2) + 3*x + 6

def function_2d(x : float, y : float) -> float:
    return 5*(x**2)*(y**2) + 3*x*y + 6

integral_1d = GaussIntegral(2, function_1d)

print(integral_1d.integrate())