import sympy as sp

# sp.init_printing()

k1 = sp.symbols('k1')
k2 = sp.symbols('k2')
k3 = sp.symbols('k3')

# v = sp.symbols('v')
# u = sp.symbols('u')
u = (1 + sp.sqrt(1 + 4*k3))/(2*k3)
v = (1 + u)/(k1 * u)

# M = sp.Matrix([
#     [
#         1 - k1*v*(1/(1+u)**2),
#         -k1*u**2/(1+u),
#     ],
#     [
#         k2*k3*v*(2*u + u**2)/(1+u)**2,
#         k2*(1 - k3*u**2/(1+u)),
#         ]
# ])
# sp.pretty_print(M)
lam = sp.symbols('lambda')
# cp = sp.det(M - lam * sp.eye(2))
# eigs = sp.roots(sp.Poly(cp, lam))
# print eigs

f1 = k3*u**2/(1+u) + k1*v/(1+u)**2 - 2

f2 = 1 - k1*v/(1+u)**2 - k3*u**2/(1+u) + k1*k3*v*u**2/(1+u)**3 + k1*k3*u**3*(2+u)/(1+u)**3
# eigs = sp.roots(sp.Poly(lam**2 + lam*f1 + f2, lam))
eigs = sp.solve(lam**2 + lam*f1 + f2, lam)
print eigs

