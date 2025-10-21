# Importar as bibliotecas necessárias
# Definir os valores simbólicos das propriedades do material
# Cálculo de v21 a
# Calcular as rigidezes de flexão equivalentes
# Definir o intervalo de valores de x e y e criar uma malha bidimensional
# Calcular a deflexão w(x, y) somando os termos das séries duplas de Fourier
# com os índices m e n ímpares (1, 3, 5)
# Representar graficamente a deflexão obtida em 3D
# Determinar os valores da deflexão máxima e as coordenadas correspondentes

import math
import numpy as np
import plotly.graph_objects as go

#Dados relativos ao caso 6
a= 0.75   #m
b= 0.6    #m
t= 0.002  #m
E1= 130e9 #Pa
E2= 10e9  #Pa
v12= 0.26
#Das aulas de estruturas 2 temos que :
v21 = (E2 / E1) * v12
G12= 5e9  #Pa
p0= -175   #N/mm^2
pi=math.pi

#Da alinea a) anteiror obtemos as seguintes rigidezes
D_x= (E1*t**3)/(12*(1-v12*v21))
D_y= (E2*t**3)/(12*(1-v12*v21))
D_xy=(G12*t**3)/12


x = np.linspace(0, a, 120)
y = np.linspace(0, b, 120)
X, Y = np.meshgrid(x, y)
w = np.zeros_like(X)


valores_m= [1,3,5]
valores_n= [1,3,5]

for m in valores_m:
    for n in valores_n:
        p_mn = (16 * p0) / (pi**2 * m * n)
        denominador = (
            D_x * ((m * pi / a)**4)
            + (D_x * v21 + 4 * D_xy + D_y * v12) * ((m * pi / a)**2 * (n * pi / b)**2)
            + D_y * ((n * pi / b)**4)
        )
        a_mn = p_mn / denominador
        w += a_mn * np.sin(m * pi * X / a) * np.sin(n * pi * Y / b) #como é um somatório += serve para guardar todo os valroes e somar-los


fig = go.Figure(data=[
    go.Surface(
        x=X, y=Y, z=w,
        colorscale='Turbo',
        contours_z=dict(show=True, usecolormap=True, project_z=True)
    )
])

fig.update_layout(
    title="Deflexão w(x, y) da placa ortotrópica",
    scene=dict(
        xaxis_title="x [m]",
        yaxis_title="y [m]",
        zaxis_title="w [m]",
        aspectratio=dict(x=1, y=1, z=0.3)
    ),
    width=900,
    height=650,
    margin=dict(l=0, r=0, t=50, b=0)
)

fig.show()

w_min = np.min(w)
indice_min = np.unravel_index(np.argmin(w), w.shape)
x_min = X[indice_min]
y_min = Y[indice_min]

print(f"Deflexão máxima (em valor absoluto): {w_min:.3e} m")
print(f"Ocorre em: x = {x_min:.3f} m, y = {y_min:.3f} m")
