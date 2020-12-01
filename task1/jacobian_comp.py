import sympy as sp

# sp.init_printing()

k1 = sp.symbols('k1')
k2 = sp.symbols('k2')
# k3 = sp.symbols('k3')
# k1 = 100
# k2 = -100
k3 = 1

# v = sp.symbols('v')
# u = sp.symbols('u')
u = (1 + sp.sqrt(1 + 4 * k3)) / (2 * k3)
v = (1 + u) / (k1 * u)

# M = sp.Matrix([
#     [
#         1 - k1*v*(1/(1+u)**2),
#         -k1*u**2/(1+u),
#     ],
#     [
#         -k2*k3*v*(2*u + u**2)/(1+u)**2,
#         k2*(1 - k3*u**2/(1+u)),
#         ]
# ])
# sp.pretty_print(M)
lam = sp.symbols('lambda')

eigs = sp.solve(
    (1 - k1 * v * (1 / (1 + u) ** 2) - lam) * (k2 * (1 - k3 * u ** 2 / (1 + u)) - lam) - (-k1 * u ** 2 / (1 + u)) * (
                -k2 * k3 * v * (2 * u + u ** 2) / (1 + u) ** 2), lam)

sp.print_latex(eigs[0])
print eigs[0].evalf()
sp.print_latex(eigs[1])
print eigs[1].evalf()

print sp.print_latex(
    (1 - k1 * v * (1 / (1 + u) ** 2) - lam) * (k2 * (1 - k3 * u ** 2 / (1 + u)) - lam) - (-k1 * u ** 2 / (1 + u)) * (
                -k2 * k3 * v * (2 * u + u ** 2) / (1 + u) ** 2))

# print sp.print_latex((-1 - k2 + k1*v/(u**2 +2*u + 1) + k2*k3*u*u/(u+1)))

# sp.print_latex(sp.expand(-4*(-k1*k2*k3*u**4*v/(u**3 + 3*u**2 + 3*u + 1) - 2*k1*k2*k3*u**3*v/(u**3 + 3*u**2 + 3*u + 1) + k1*k2*k3*u**2*v/(u**3 + 3*u**2 + 3*u + 1) - k1*k2*v/(u**2 + 2*u + 1 )- k2*k3*u**2/(u + 1)+ k2) + (k1*v/(u**2 + 2*u + 1) + k2*k3*u**2/(u + 1)  - k2 - 1)**2))

print (0.381966011250105/1.17557050458481)**2
      # - 0.105572809000084