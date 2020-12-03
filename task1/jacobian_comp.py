import sympy as sp


k1 = sp.symbols('k1')
k2 = sp.symbols('k2')
# k3 = sp.symbols('k3')
k3 = 1

# v = sp.symbols('v')
# u = sp.symbols('u')
u = (1 + sp.sqrt(1 + 4 * k3)) / (2 * k3)
v = (1 + u) / (k1 * u)

lam = sp.symbols('lambda')

#2-d system

# eigs = sp.solve(
#     (1 - k1 * u * v * ((2+u) / ((1 + u) ** 2)) - lam) * (k2 * (1 - k3 * u ** 2 / (1 + u)) - lam) - (-k1 * u ** 2 / (1 + u)) * (
#                 -k2 * k3 * v * (2 * u + u ** 2) / (1 + u) ** 2), lam)
#
# sp.print_latex(eigs[0])
# print eigs[0].evalf()
# sp.print_latex(eigs[1])
# print eigs[1].evalf()

# print sp.print_latex(
#     (1 - k1 * u * v * (2+u) / ((1 + u) ** 2) - lam) * (k2 * (1 - k3 * u ** 2 / (1 + u)) - lam) - (-k1 * u ** 2 / (1 + u)) * (
#               -k2 * k3 * v * (2 * u + u ** 2) / (1 + u) ** 2))
# print (0.190983005625053 /1.17557050458504)**2 - 0.026393202250021


#3-d system

c = sp.symbols('c')

J11 = (1 - k1 * u * v * (2+u) / ((1 + u) ** 2))/c
J12 = (-k1 * u ** 2 / (1 + u))/c
J31 = k2 * k3 * v * (2 * u + u ** 2) / (1 + u) ** 2
J32 = (-k2 * (1 - k3 * u ** 2 / (1 + u)))

# print sp.print_latex(J11)
# print sp.print_latex(J12)
# print sp.print_latex(J31)
# print sp.print_latex(J32)

# characteristic_polynomial = (J11 - lam)*(-lam)*(c-lam) + J31*J12*(c-lam) + 0 - (J11 - lam)*J32
#
# print sp.print_latex(characteristic_polynomial)
# eigs_2 = sp.solve(characteristic_polynomial, lam)
# sp.print_latex(eigs_2[0])
# print (sp.expand(eigs_2[0])).evalf()
# sp.print_latex(eigs_2[1])
# print (sp.expand(eigs_2[1])).evalf()
# sp.print_latex(eigs_2[2])
# print (sp.expand(eigs_2[2])).evalf()


characteristic_polynomial_o = (1/c - lam)*(-lam)*(c-lam) - (1/c - lam)*(-k2)
print sp.print_latex(characteristic_polynomial_o)
eigs_2_o = sp.solve(characteristic_polynomial_o, lam)
sp.print_latex(eigs_2_o[0])
print (sp.expand(eigs_2_o[0])).evalf()
sp.print_latex(eigs_2_o[1])
print (sp.expand(eigs_2_o[1])).evalf()
sp.print_latex(eigs_2_o[2])
print (sp.expand(eigs_2_o[2])).evalf()