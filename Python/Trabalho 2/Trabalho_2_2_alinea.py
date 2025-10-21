# Importar bibliotecas
# Devido à complexidade do código, vamos dividir o mesmo por etapas
# --------------------------
# - Definição simbólica das constantes, coordenadas e dos parâmetros físicos
# - Construção da solução particular e da solução homogênea da deflexão w(x,y)
# - Definição total da equação da deflexão -> w(x,y)
# --------------------------
# - Derivadas parciais em x e em y 
# --------------------------
# - Cálculo dos momentos e dos esforços de corte nos limites y=0 e y=b
# --------------------------
# - Definição das condições de fronteira
# - Simplificação dividindo pelo fator sin(alpha*x)
# --------------------------
# - Preparação para obter a solução simbólica de Am, Bm, Cm e Dm
# --------------------------
# - Definição das propriedades da placa e parâmetros do problema
# - Cálculo de Pm para cada modo m
# - Substituição numérica e resolução do sistema linear para cada modo
# -----------------------------
# - Cálculo da deflexão total w(x,y) somando todos os modos m
# -----------------------------
# - Criar o gráfico 3D 
# -----------------------------
# - Determinar os valores de w no centro da placa e as respetivas coordenadas
# -----------------------------

import sympy as sp
import numpy as np
import plotly.graph_objects as go

# -------------------------------------------------------------
# 1) Símbolos e solução simbólica
# -------------------------------------------------------------
Am, Bm, Cm, Dm = sp.symbols('Am Bm Cm Dm')
x, y, a, b, alpha, nu, D, Pm, kw, ktheta = sp.symbols('x y a b alpha nu D Pm kw ktheta', real=True)

W_part = (Pm / D) / alpha**4
W_h = Am*sp.cosh(alpha*y) + Bm*sp.sinh(alpha*y) + Cm*y*sp.cosh(alpha*y) + Dm*y*sp.sinh(alpha*y)

W_total = W_h + W_part

w_xy = W_total * sp.sin(alpha*x)

w_y   = sp.diff(w_xy, y)
w_yy  = sp.diff(w_y, y)
w_yyy = sp.diff(w_yy, y)
w_xx  = sp.diff(w_xy, x, 2)

subs0 = {y: 0}
subsb = {y: b}

w0, w_y0, w_yy0, w_yyy0 = [sp.simplify(expr.subs(subs0)) for expr in [w_xy, w_y, w_yy, w_yyy]]
w_b, w_yb, w_yyb, w_yyyb = [sp.simplify(expr.subs(subsb)) for expr in [w_xy, w_y, w_yy, w_yyy]]

M_y0 = D * (w_yy0 - nu * w_xx.subs(subs0))
M_yb = D * (w_yyb - nu * w_xx.subs(subsb))
Q_y0 = D * (w_yyy0 - (2 - nu) * alpha**2 * w_y0)
Q_yb = D * (w_yyyb - (2 - nu) * alpha**2 * w_yb)

eq1_x = sp.simplify(Q_y0 + kw*w0)
eq2_x = sp.simplify(M_y0)
eq3_x = sp.simplify(Q_yb - kw*w_b)
eq4_x = sp.simplify(M_yb - ktheta*w_yb)

sin_factor = sp.sin(alpha*x)
eqs = [sp.simplify(eq / sin_factor) for eq in [eq1_x, eq2_x, eq3_x, eq4_x]]

unknowns = [Am, Bm, Cm, Dm]
M_sym = sp.zeros(4,4)
RHS_sym = sp.zeros(4,1)

for i, eq in enumerate(eqs):
    for j, u in enumerate(unknowns):
        M_sym[i,j] = sp.diff(eq, u)
    const_term = eq.subs({u:0 for u in unknowns})
    RHS_sym[i] = -const_term

sol = sp.solve(M_sym * sp.Matrix(unknowns) - RHS_sym, unknowns, dict=True)[0]
A_expr = sp.simplify(sol[Am])
B_expr = sp.simplify(sol[Bm])
C_expr = sp.simplify(sol[Cm])
D_expr = sp.simplify(sol[Dm])

print('--- Soluções simbólicas gerais (simplificadas) ---')
print(f"A_m = {A_expr}")
print(f"B_m = {B_expr}")
print(f"C_m = {C_expr}")
print(f"D_m = {D_expr}\n")

a_val = 0.75
b_val = 0.6
t_val = 0.001
E_val = 210e9
nu_val = 0.30
kw_val = 125e3
ktheta_val = 1060
p0_val = 175
D_val = E_val * t_val**3 / (12 * (1 - nu_val**2))

def compute_Pm_alpha(m, p0, a):
    return (4.0 * p0 / (a * (m * np.pi / a))), (m * np.pi / a)

m_list = [1, 3, 5, 7, 9]
ABCD = {}

for m in m_list:
    Pm_val, alpha_val = compute_Pm_alpha(m, p0_val, a_val)
    subs_vals = {
        alpha: alpha_val, a: a_val, b: b_val, D: D_val, nu: nu_val,
        kw: kw_val, ktheta: ktheta_val, Pm: Pm_val
    }

    M_num = np.array(sp.N(M_sym.subs(subs_vals)), dtype=float)
    RHS_num = np.array(sp.N(RHS_sym.subs(subs_vals)), dtype=float).flatten()

    try:
        sol_num = np.linalg.solve(M_num, RHS_num)
    except np.linalg.LinAlgError:
        sol_num = np.linalg.pinv(M_num) @ RHS_num

    ABCD[m] = dict(A=sol_num[0], B=sol_num[1], C=sol_num[2], D=sol_num[3],
                   Pm=Pm_val, alpha=alpha_val)

nx, ny = 100, 100
x_vals = np.linspace(0, a_val, nx)
y_vals = np.linspace(0, b_val, ny)
X, Y = np.meshgrid(x_vals, y_vals)
W = np.zeros_like(X)

for m in m_list:
    s = ABCD[m]
    alpha_m, Pm_m = s['alpha'], s['Pm']
    A, B, C, Dm_ = s['A'], s['B'], s['C'], s['D']
    cosh_term = np.cosh(alpha_m * Y)
    sinh_term = np.sinh(alpha_m * Y)
    Wm = (A * cosh_term + B * sinh_term + C * Y * cosh_term + Dm_ * Y * sinh_term + Pm_m / (D_val * alpha_m**4))
    W += Wm * np.sin(alpha_m * X)

W = -W

fig = go.Figure(data=[
    go.Surface(
        z=W,
        x=X,
        y=Y,
        colorscale='Turbo',
        colorbar=dict(title='w (m)'),
        contours={"z": {"show": True, "usecolormap": True, "project_z": True}}
    )
])

fig.update_layout(
    title='Deflexão w(x,y) - Caso 6 (Lévy)',
    scene=dict(
        xaxis_title='x (m)',
        yaxis_title='y (m)',
        zaxis_title='w (m)'
    ),
    width=900,
    height=700
)

fig.show()

ix = np.argmin(abs(x_vals - a_val/2))
iy = np.argmin(abs(y_vals - b_val/2))
w_centro = W[iy, ix]
print(f"w_centro ≈ {w_centro:.3e} m (negativo = para baixo)")

# -------------------------------------------------------------
# Passar os valores obtidos para um documento de texto 
# -------------------------------------------------------------
fname = 'resultados_Levy_Caso6.txt'
with open(fname, 'w') as f:
    f.write('Resultados - Plate Bending (Lévy) - Caso 6\n')
    f.write('------------------------------------------\n')
    f.write(f'a = {a_val:.3f} m, b = {b_val:.3f} m, t = {t_val:.4f} m\n')
    f.write(f'E = {E_val:.2e} Pa, nu = {nu_val:.2f}\n')
    f.write(f'kw = {kw_val:.2e} N/m^3, ktheta = {ktheta_val:.2f} N·m/rad, p0 = {p0_val:.2f} N/m^2\n')
    f.write(f'D = {D_val:.3e} N·m\n\n')
    f.write('--- Soluções simbólicas gerais ---\n')
    f.write(f'A_m = {A_expr}\nB_m = {B_expr}\nC_m = {C_expr}\nD_m = {D_expr}\n\n')
    f.write('Constantes por modo:\n')
    f.write(' m |       A [m]        |       B [m]        |       C [m]        |       D [m]\n')
    f.write('---|--------------------|--------------------|--------------------|--------------------\n')
    for m in m_list:
        s = ABCD[m]
        f.write(f'{m:2d} | {s["A"]:+.5e} | {s["B"]:+.5e} | {s["C"]:+.5e} | {s["D"]:+.5e}\n')
    f.write(f'\nDeflexão no centro:\n w_centro = {w_centro:.5e} m (negativo = para baixo)\n')

print(f'Resultados exportados para "{fname}" com sucesso.')
