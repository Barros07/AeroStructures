#Defenir todas as constantes que temos 
#Resolver a equação diferencial EI*v''''(x) = -q
#Integrar 4 vezes e obter c1, c2, c3 e c4
#Defenir as condições de fronteira
#Substituir e resolver o sistema para calcular os valores das constantes c1, c2, c3 e c4 através das condições de fronteira
#Obter a equação da linha elástica v(x)
#Obter a equação do momento fletor no centro da viga (x=L/2)
#Obter a equação da força vertical no ponto mais a direita (x=L)
#Substituir pelos valores numerios as constantes defenidas na linha 1
#Calcular a deflexão (ja temos linha 6) 
#calcular a inclinação= theta(x)= v'(x)
#Desenhar graficamente os dados da linha 10 e 11 ao longo de 0<x<L

#Bibliotecas necessárias para a resoluçao do exercicio
import sympy as sp #calculo diferencial
import numpy as np #transformar os dados e analisar todos os pontos para ser possível desenhar os gráficos
import matplotlib.pyplot as plt #Representar gráficamente as o declive e a deflexão ao longo da viga

#Defenir as variáveis
x = sp.symbols ('x')
E, I, L, q, k=sp.symbols('E I L q k')
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
#Defenir a equaçao diferencial
v= sp.Function('v')(x)
EqDif= sp.Eq(E*I*sp.diff(v,(x,4)),-q)

#Resolver a EqDif
v_solv= sp.dsolve(EqDif)
V_geral= v_solv.rhs #Obter a expressão da linha elástica
print ("\n Equação da linha elástica (com constantes)")
sp.pprint(V_geral)

#Vamos defenir as condições de fronteira (4 equaçoes de fronteira, 2 para cada extermidade)
#Devido ao encastramento
 #v(x)=0 --- vai anular c4- condição 1
 #v'(x)=0---- vai anular c3- condição 2
#Devido ao apoio da mola
#v(L)=0--- vamos obter uma equação dependente de c1 e c2-condição 3
#Equação do momento
#Para a eq4 -> M_L+M_mola=0 <=> M_L= -k*v'(L) <=> EI*v''(L)= K*v'(L)<=> EI*v''(L)-K*v'(L)=0- condição 4 
cond1 = sp.Eq(V_geral.subs(x,0), 0)
cond2 = sp.Eq(sp.diff(V_geral,x).subs(x,0), 0)
cond3 = sp.Eq(V_geral.subs(x,L), 0)
cond4 = sp.Eq(E*I*sp.diff(V_geral,(x,2)).subs(x,L) - k*sp.diff(V_geral,x).subs(x,L), 0)

#Resolver sistema para c1,c2,c3 e c4
sol_const = sp.solve([cond1, cond2, cond3, cond4], (C1, C2, C3, C4))
print("\nConstantes de integração encontradas:")
print(sol_const)

#Substituir na equação final da linha elástica
v_final = V_geral.subs(sol_const)

print("\nEquação da linha elástica v(x):")
sp.pprint(sp.simplify(v_final))

#Inclinação, momento e força vertical
theta= sp.diff(v_final,x)
M= E*I*sp.diff(v_final,(x,2)) 
F= -E*I*sp.diff (v_final,(x,3))  #equivalente à 3ª derivada da deflexao

print("\n Inclinaçao :")
sp.pprint(sp.simplify(theta))
print("\n Momento:")
sp.pprint(sp.simplify(M))
print("\n Força vertical:")
sp.pprint(sp.simplify(F))

#Introduzir os valores do caso 6
b_val= 0.01
h_val= 0.01
L_val= 0.8
E_val= 70e9
q_val= 250
k_val= 300
I_val= (b_val * h_val**3)/12

#Substituir os valores nas equações que tivemos a defenir até agora 
subs_dict = {E: E_val, I: I_val, L: L_val, q: q_val, k: k_val}
v_val, theta_val, M_val, F_val = [
    sp.simplify(expr.subs(subs_dict)) 
    for expr in (v_final, theta, M, F)
]

print("\n A linha elastica para o caso ......, v(x)")
sp.pprint(v_val)
print("\n A inclinaçao para o caso ......, theta(x)")
sp.pprint(theta_val)
print("\n A linha elastica para o caso ......, M(x)")
sp.pprint(M_val)
print("\n A linha elastica para o caso ......, F(x)")
sp.pprint(F_val)

#Transformar os valores de função simbólica para uma função numpy
def evaluate_expr(expr, x_vals):

    f = sp.lambdify(x, expr, "numpy")
    try:
        y = f(x_vals) 
        y_arr = np.asarray(y)
        if y_arr.dtype == object:
            y_arr = np.array([float(sp.N(val)) for val in y_arr], dtype=float)
        else:
            y_arr = y_arr.astype(float)
        return y_arr
    except Exception:
        return np.array([float(sp.N(expr.subs(x, float(t)))) for t in x_vals], dtype=float)

#De acordo com o caso que temos ,vamos calcular o momento a meio da viga e a força de apoio em L
M_meio = float(sp.N(M_val.subs(x, L_val/2)))
F_dir = float(sp.N(F_val.subs(x, L_val)))

#Apresentar os resultados da alínea a)

print("\nOs resultados da alínea a) são")
print(f"Momento a meio da viga M(L/2) [N·m]: {M_meio:.2f}")
print(f"Força no apoio direito F(L) [N]: {F_dir:.4f}")

#Representaçao da alínea b)
subs_dict = {E: E_val, I: I_val, L: L_val, q: q_val, k: k_val}
v_num = sp.simplify(v_final.subs(subs_dict))
theta_num = sp.simplify(theta.subs(subs_dict))

# Criar funções numéricas com lambdify
v_fun = sp.lambdify(x, v_num, "numpy")
theta_fun = sp.lambdify(x, theta_num, "numpy")

#Incremento dos valores (maior numero=maior precisão)
x_valores = np.linspace(0, float(L_val), 100)

#Substituir pelos valores de x da função
v_valores = evaluate_expr(v_num, x_valores)
theta_valores = evaluate_expr(theta_num, x_valores)

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
plt.plot(x_valores, v_valores, label="Deflexão (v)", color="orange")
plt.title("Deflexão ao longo da viga")
plt.xlabel("Comprimento da viga (m)")
plt.ylabel("Deflexão (m)")
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(x_valores, theta_valores, label="Inclinação (v')", color="blue")
plt.title("Inclinação ao longo da viga")
plt.xlabel("Comprimento da viga (m)")
plt.ylabel("Inclinação (rad)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

caminho_saida = r"C:\Users\HP\Desktop\ubi\2025-2026\1º Semestre\Placas e Cascas\Trabalho-1\resultados_viga.txt"
with open(caminho_saida, "w", encoding="utf-8") as f:
    f.write("===== RESULTADOS DO CÁLCULO DA VIGA =====\n\n")

    f.write("1) Equação da linha elástica geral (com constantes):\n")
    f.write(str(V_geral) + "\n\n")

    # <<< NOVO BLOCO PARA GUARDAR AS CONSTANTES DE INTEGRAÇÃO >>>
    f.write("2) Constantes de integração (C1..C4):\n")
    for const, val in sol_const.items():
        f.write(f"{const} = {sp.simplify(val)}\n")
    f.write("\n")

    f.write("3) Equação da linha elástica v(x) (com C1..C4 substituídos):\n")
    f.write(str(sp.simplify(v_final)) + "\n\n")

    f.write("4) Inclinação θ(x):\n")
    f.write(str(sp.simplify(theta)) + "\n\n")

    f.write("5) Momento M(x):\n")
    f.write(str(sp.simplify(M)) + "\n\n")

    f.write("6) Força cortante F(x):\n")
    f.write(str(sp.simplify(F)) + "\n\n")

    f.write("===== RESULTADOS COM VALORES NUMÉRICOS =====\n\n")
    f.write(f"Dimensões: b = {b_val} m, h = {h_val} m, L = {L_val} m\n")
    f.write(f"Constantes: E = {E_val}, q = {q_val}, k = {k_val}, I = {I_val}\n\n")

    f.write("v(x) numérico:\n")
    f.write(str(v_val) + "\n\n")

    f.write("θ(x) numérico:\n")
    f.write(str(theta_val) + "\n\n")

    f.write("M(x) numérico:\n")
    f.write(str(M_val) + "\n\n")

    f.write("F(x) numérico:\n")
    f.write(str(F_val) + "\n\n")

    f.write("===== VALORES ESPECIAIS =====\n\n")
    f.write(f"Momento a meio da viga M(L/2) [N*m]: {M_meio:.2f}\n")
    f.write(f"Força no apoio direito F(L) [N]: {F_dir:.5f}\n")
