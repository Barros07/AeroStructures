# ==================================
# MÉTODO DE RITZ — PLACA ORTOTRÓPICA
# ==================================
#Introduçao de bibliotecas
#Introdução dos parâmetros geométricos 
#Calculo das rigidezes de flexão
#Pedir o numero de modos
#Definição da função X e Y
#Definiçao da funçao do deslocamento w
#Definição das matrizes R e B
#Resolver o problema generalizado de autovalores
#Pedir ao utilzador o modo a visualizar
#Reconstruir a funçao w(x,y)
#Reprensentação simbólica da frequência omega
#Cálculo do deslocamento máximo |w|max
#Representação gráfica do modo de vibração

import sympy as sp                   # Biblioteca simbólica (derivadas e integrais)
import numpy as np                   # Biblioteca numérica (matrizes, vetores)
from scipy.linalg import eigh        # Resolver problemas generalizados de autovalores
import matplotlib.pyplot as plt    # Biblioteca para gráficos 3D interativos


a = 0.750     # [m]
b = 0.600     # [m]
t = 0.002     # [m]
rho = 1600    # [kg/m³]
E1 = 130e9    # [Pa]
E2 = 10e9     # [Pa]
nu12 = 0.26
G12 = 5e9     # [Pa]

nu21 = nu12 * E2 / E1

Q11 = E1 / (1 - nu12 * nu21)
Q22 = E2 / (1 - nu12 * nu21)
Q12 = nu12 * E2 / (1 - nu12 * nu21)
Q66 = G12

D11 = Q11 * t**3 / 12
D22 = Q22 * t**3 / 12
D12 = Q12 * t**3 / 12
D66 = Q66 * t**3 / 12

try:
    n_modos = int(input("Introduz o número de modos de vibração: "))
except ValueError:
    n_modos = 2

if n_modos < 1:
    n_modos = 2

M = n_modos
N = n_modos
n = M * N   # número total de graus de liberdade



# ============================
x, y = sp.symbols('x y', real=True)
X = []
Y = []

for i in range(1, M + 1):
    X_i = (x / a)**(i + 1) - 2 * (x / a)**(i + 2) + (x / a)**(i + 3)
    X.append(X_i)

for j in range(1, N + 1):
    Y_j = (y / b)**j - (y / b)**(j + 1)
    Y.append(Y_j)

c = sp.symbols(f'c0:{n}', real=True)
w = 0
k = 0

for i in range(M):
    for j in range(N):
        w += c[k] * X[i] * Y[j]
        k += 1


R = np.zeros((n, n), dtype=float)
B = np.zeros((n, n), dtype=float)

for p in range(n):
    wp = sp.diff(w, c[p])
    wp_xx = sp.diff(wp, x, 2)
    wp_yy = sp.diff(wp, y, 2)
    wp_xy = sp.diff(sp.diff(wp, x), y)

    for q in range(n):
        wq = sp.diff(w, c[q])
        wq_xx = sp.diff(wq, x, 2)
        wq_yy = sp.diff(wq, y, 2)
        wq_xy = sp.diff(sp.diff(wq, x), y)

        expr_R = (D11 * wp_xx * wq_xx +
                  2 * D12 * wp_xx * wq_yy +
                  D22 * wp_yy * wq_yy +
                  4 * D66 * wp_xy * wq_xy)

        expr_B = rho * t * wp * wq

        R[p, q] = float(sp.integrate(sp.integrate(expr_R, (x, 0, a)), (y, 0, b)))
        B[p, q] = float(sp.integrate(sp.integrate(expr_B, (x, 0, a)), (y, 0, b)))


omega2, vec = eigh(R, B)
omega2 = np.real(omega2[omega2 > 1e-8])
omega = np.sqrt(omega2)
freq = omega / (2 * np.pi)

print("\n=== FREQUÊNCIAS NATURAIS ===")
for i, f in enumerate(freq, start=1):
    print(f"Modo {i}: f = {f:.3f} Hz")


try:
    modo = int(input("\nIntroduz o número do modo a visualizar (ex: 1): "))
except ValueError:
    modo = 1

if modo < 1 or modo > len(freq):
    modo = 1

phi = vec[:, modo - 1]
phi = phi / np.max(np.abs(phi))


Nx, Ny = 80, 80
x_vals = np.linspace(0, a, Nx)
y_vals = np.linspace(0, b, Ny)
Xgrid, Ygrid = np.meshgrid(x_vals, y_vals)
W = np.zeros_like(Xgrid, dtype=float)

k = 0
for i in range(M):
    for j in range(N):
        X_fun = sp.lambdify(x, X[i], "numpy")
        Y_fun = sp.lambdify(y, Y[j], "numpy")
        W += phi[k] * X_fun(Xgrid) * Y_fun(Ygrid)
        k += 1


x, y, c1, a_s, b_s, t_s, rho_s, D11_s, D22_s, D12_s, D66_s = sp.symbols(
    'x y c1 a b t rho D11 D22 D12 D66', real=True)

X1 = (x / a_s)**2 - 2 * (x / a_s)**3 + (x / a_s)**4
Y1 = (y / b_s) - (y / b_s)**2
w_sym = c1 * X1 * Y1


w_xx = sp.diff(w_sym, x, 2)
w_yy = sp.diff(w_sym, y, 2)
w_xy = sp.diff(sp.diff(w_sym, x), y)

U = sp.integrate(sp.integrate(D11_s*w_xx**2 + 2*D12_s*w_xx*w_yy +
                              D22_s*w_yy**2 + 4*D66_s*w_xy**2, (y, 0, b_s)), (x, 0, a_s))
T = sp.integrate(sp.integrate(rho_s * t_s * w_sym**2, (y, 0, b_s)), (x, 0, a_s))

omega_symbolic = sp.simplify(sp.sqrt(U / T))

print("\n=== EXPRESSÃO SIMBÓLICA DA FREQUÊNCIA ω ===")
sp.pretty_print(omega_symbolic)
print("\n=== EXPRESSÃO EM LATEX ===")
print(sp.latex(omega_symbolic))

w_max = np.max(np.abs(W))
print(f"\nValor máximo do deslocamento relativo |w|max = {w_max:.4f}")

fig = plt.figure(facecolor='white')  # fundo da figura branco
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('white')  # fundo do gráfico branco

# Superfícielate
ax.plot_surface(Xgrid, Ygrid, -W, cmap='turbo', edgecolor='none')

# Rótulos dos eixos
ax.set_xlabel('x (m)', color='black')
ax.set_ylabel('y (m)', color='black')
ax.set_zlabel('w(x,y)', color='black')

# Título
ax.set_title(f"Modo {modo} — f = {freq[modo-1]:.3f} Hz  |w|max = {w_max:.4f}", color='black')

# Cores dos rótulos e ticks
ax.xaxis.label.set_color('black')
ax.yaxis.label.set_color('black')
ax.zaxis.label.set_color('black')
ax.tick_params(colors='black')

plt.show()
