# Importar as bibliotecas
# Definir as constantes simbólicas que representam as propriedades do material
# Definir a deflexão em função de x e y → w(x, y)
# Definir os coeficientes de rigidez reduzida
# Definir a matriz de rigidez reduzida Q para um material ortotrópico
# Determinar as curvaturas da superfície média da placa
# Definir as tensões em função das deformações
# Cálculo simbólico integral das tensões para obter os momentos
# Simplificando o resultado obtido, obtemos os valores das constantes de rigidez de uma placa ortotrópica

import sympy as sp
import math 
import numpy as np

E1,E2,v12,v21,G12=sp.symbols('E1,E2,v12,v21,G12')
sigma_x, sigma_y, tau_xy= sp.symbols('sigma_x, sigma_y, tau_xy')
ep_x,ep_y,gamma_xy=sp.symbols('ep_x,ep_y,gamma_xy')
x,y,z,t= sp.symbols('x,y,z,t')
w=sp.Function('w')(x,y)

Q11=E1/(1-v12*v21)
Q22=(E2/E1)*Q11
Q12=v12*Q11
Q66= G12

Q=sp.Matrix([
    [Q11, Q12, 0],
    [Q12,Q22,0],
    [0,0,Q66]
])

k_xx=sp.diff(w,x,2)
k_yy=sp.diff(w,y,2)
k_xy= 2*sp.diff(w,x,y)

k=sp.Matrix([k_xx,k_yy,k_xy])

eps=-z*k

tensoes= Q*eps

M = sp.integrate(z * tensoes, (z, -t/2, t/2))

D = sp.integrate(z**2 * Q, (z, -t/2, t/2))
D_simplificado = sp.simplify(D)

M_relaçao = sp.Eq(M, D_simplificado * k)

sp.pprint(D_simplificado)