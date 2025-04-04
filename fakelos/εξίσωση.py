#%% Λύση δευτεροβάθμιας εξίσωσης 0.04x^2 - 0.1184x - 0.96 = 0

from sympy import symbols, Eq, solve

# Ορισμός μεταβλητής και εξίσωσης
x = symbols('x')
quadratic_eq = Eq(0.04 * x**2 - 0.1184 * x - 0.96, 0)

# Επίλυση εξίσωσης
quadratic_solutions = solve(quadratic_eq, x)
quadratic_solutions

# %%
