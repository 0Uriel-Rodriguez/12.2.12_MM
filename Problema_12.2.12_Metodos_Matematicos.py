#!/usr/bin/env python
# coding: utf-8

# In[62]:


import numpy as np
from scipy.special import legendre

n = 2  # Cambia n según lo que necesites
coeff = legendre(n)  # Calcula los coeficientes de Legendre para n
roots = np.roots(coeff)  # Encuentra las raíces

# Filtra las raíces reales positivas
positive_real_roots = [root for root in roots if root.real >= 0 and np.iscomplex(root) == False]

# Encuentra la raíz real positiva más grande (mayor valor)
if positive_real_roots:
    largest_positive_real_root = max(positive_real_roots)
    print(f"Raíz real positiva más grande para n = {n}:")
    print(largest_positive_real_root.real)
else:
    print(f"No se encontraron raíces positivas reales para n = {n}.")


# In[61]:


import sympy as sp

# Símbolo x como variable
x = sp.symbols('x')

# Grado del polinomio de Legendre
n = 2  # Cambia n según lo que necesites

# Expresión hipergeométrica de Pn(x)
hypergeometric_expression = sp.hyper((-n, n + 1), (1,), 1/2 * (1 - x))

# Definir la función F(x) utilizando la representación hipergeométrica
F = (x**2 - 1)**n * hypergeometric_expression

# Derivada de F(x) con respecto a x
F_prime = sp.diff(F, x)

# Método de Newton-Raphson para encontrar la raíz más grande
x0 = 0.5  # Estimación inicial cercana a la raíz
tolerance = 1e-8  # Tolerancia para convergencia
max_iterations = 100  # Número máximo de iteraciones

for i in range(max_iterations):
    # Evaluación de F(x) y su derivada en x0
    f_x0 = F.evalf(subs={x: x0})
    f_prime_x0 = F_prime.evalf(subs={x: x0})

    # Actualización de x0 utilizando el método de Newton-Raphson
    x0 = x0 - f_x0 / f_prime_x0

    # Verificación de convergencia
    if abs(f_x0) < tolerance:
        break

# La raíz más grande se encuentra en x0
print(f"La aproximación de la raíz más grande de P_{n}(x) es: {x0:.10f}")


# In[ ]:





# In[ ]:




