# ==================================
# MÉTODO C/H SIMBÓLICO + CONVERGÊNCIA
#===================================
#Importar bibliotecas
#Definir as porpriedades do material
#Definir a matriz D para o caso ortotrópico
#Definir o tipo de borda em cada lado
#Definição simbolica das variaveis
    #Definir a matriz H simbolica e inicializaçao em 0
    #Definir o polinomio w_e (a1...a12)
    #Definir as rotações theta_x e theta_y (derivadas)
    #construir a matriz simbolica C e inicializaçao em 0
    #Preencher a matriz C simbolica com ([1,4,7,10] -> w; [2,5,8,11] -> theta_x; [3,6,9,12] -> theta_y)
    #Integrar simbolicamente I_xy = ∫∫ H' * D * H dx dy
#Iniciar o metodo de convergencia
    #Definir o erro relativo de convergencia
    #Iniciar a malha com 1 quadrado por lado
    #iniciar o loop de convergencia
        #Incrementar o numero de quadrados por lado
        #Gerar a malha de nos (cada no tem 3 graus de liberdade) e a conectividade
        #Montar a matriz K e o vetor F e inicializa-los em 0
        # Loop pelos elementos para calcular Ke e Fe e espalhar para K, F
        # Obter coordenadas do elemento
        # Avaliar integral simbolica I_xy para o elemento
        # Construir c_num
        # Substituir em cada no e preencher as linhas correspondentes
        # Verificar o condicionamento de C_num e calcular a inversa
        # Calcular a matriz elementar Ke e o vetor de cargas elementar Fe    
        #Aplicar as condições de fronteira
        #Definir a funçao para aplicar as condiçoes de fronteira (seja que tipo de fronteira for)
        #Resolver o sistema de equações lineares 
        #Montar o vetor dos deslocamentos nodais e inicializa em 0 
        #Calcular os deslocamentos nodais w
        #Calcular o erro relativo face ao passo anterior (verificar a convergencia)
#Usar a solução armnazenada
#Calcular as tensões nodais
#Resultados numéricos e impressões finais
#Represenentação gráfica dos resultados


import warnings                                       # para emitir avisos controlados
import numpy as np                                    # operações numéricas
import sympy as sp                                    # cálculo simbólico
import scipy.sparse as sp_sparse                      # matrizes esparsas
import scipy.sparse.linalg as spla                    # resolução de sistemas esparsos
from scipy.interpolate import griddata                # interpolação para grelha (griddata)
import matplotlib.pyplot as plt                       # visualização / gráficos

# ----------------------------
a = 0.625           # [m]
b = 0.650           # [m]
t = 0.001           # [m]
E1 = 135e9          # [Pa]
E2 = 10e9           # [Pa]
nu12 = 0.27         
G12 = 4e9           # [Pa]
p0 = -250.0         # [N/m^2] (negativa = para baixo)
# ----------------------------
nu21 = (E2 / E1) * nu12                             
Q11 = E1 / (1 - nu12 * nu21)                       
Q22 = E2 / (1 - nu12 * nu21)                        
Q12 = nu21 * Q11                                     
Q66 = G12                                            
Q_mat = np.array([[Q11, Q12, 0.0],
                  [Q12, Q22, 0.0],
                  [0.0, 0.0, Q66]])                 
D_b = (t**3) / 12.0 * Q_mat                  
print("\nMatriz D_b (flexão) [Pa·m³]:")
print(D_b)
# ----------------------------
B1 = 'C'   # borda x=0. encastrada
B2 = 'C'   # borda x=a. encastrada
B3 = 'C'   # borda y=0. livre
B4 = 'C'   # borda y=b. livre
# ----------------------------
x, y, x1, x2, y1, y2 = sp.symbols('x y x1 x2 y1 y2', real=True)  
D_sym = sp.Matrix(D_b)                                          
# ----------------------------
H = sp.Matrix(sp.zeros(3, 12))                                  
H[0, :] = sp.Matrix([0, 0, 0, -2, 0, 0, -6*x, -2*y, 0, 0, -6*x*y, 0]).T
H[1, :] = sp.Matrix([0, 0, 0, 0, 0, -2, 0, 0, -2*x, -6*y, 0, -6*x*y]).T
H[2, :] = sp.Matrix([0, 0, 0, 0, -2, 0, 0, -4*x, -4*y, 0, -6*x*2, -6*y*2]).T
#----------------------------
a_syms = sp.symbols('a1:13', real=True)
w_e = (
    a_syms[0]
    + a_syms[1]*x
    + a_syms[2]*y
    + a_syms[3]*x**2
    + a_syms[4]*x*y
    + a_syms[5]*y**2
    + a_syms[6]*x**3
    + a_syms[7]*x**2*y
    + a_syms[8]*x*y**2
    + a_syms[9]*y**3
    + a_syms[10]*x**3*y
    + a_syms[11]*x*y**3
)
#----------------------------
theta_x = sp.diff(w_e, x)                                       # theta_x = ∂w/∂x
theta_y = sp.diff(w_e, y)                                       # theta_y = ∂w/∂y
#----------------------------
C_sym = sp.Matrix(sp.zeros(12, 12))                          
#----------------------------
for row in range(12):
    for col in range(12):
        if row in [0, 3, 6, 9]:                                
            C_sym[row, col] = sp.diff(w_e, a_syms[col])       
        elif row in [1, 4, 7, 10]:                              
            C_sym[row, col] = sp.diff(theta_x, a_syms[col])
        else:
            C_sym[row, col] = sp.diff(theta_y, a_syms[col])
#----------------------------
print('A integrar simbolicamente H\'*D*H (apenas uma vez) ...')
Integrand = H.T * D_sym * H                                      
I_xy = sp.integrate(sp.integrate(Integrand, (x, x1, x2)), (y, y1, y2)) 
I_xy = sp.simplify(I_xy)                                        # simplificar expressão simbólica resultante
# ----------------------------
tol_rel_pct = 0.5                                               # tolerância relativa em percentagem (0.5%)
w_prev = np.nan                                                 # valor anterior w_max (para comparar)
converged = False                                               # flag de convergência
n_side = 1                                                      # inicia com 1 quadrado por lado (malha 1x1)
print(f'\nIniciando sequência de refinamento automática (critério: erro <= {tol_rel_pct:.2f}%)\n')
while not converged:
    n_side += 1                                                # incrementa número de quadrados por lado
    n_quadrados = n_side * n_side                              # número total de quadrados
    nnx = n_side + 1                                           # nós em x
    nny = n_side + 1                                           # nós em y
    n_nos = nnx * nny                                          # número total de nós
    dx = a / n_side                                            # Incremento em x
    dy = b / n_side                                            # Incremento em y

    coords = np.zeros((n_nos, 2))                              # inicializa array de coordenadas
    idx = 0                                                    # contador de nós
    for j in range(0, n_side + 1):                             # varre linhas y
        for i in range(0, n_side + 1):                         # varre colunas x
            coords[idx, :] = [i * dx, j * dy]                  # posição do nó (x, y)
            idx += 1                                           # incremento do índice

    conn = np.zeros((n_quadrados, 4), dtype=int)               # inicializa a conectividade
    ecount = 0                                                 # contador de elementos
    for j in range(1, n_side + 1):                             # varre elementos por linha
        for i in range(1, n_side + 1):                         # varre elementos por coluna
            ecount += 1                                        # índice do elemento
            n1 = (j - 1) * (n_side + 1) + i                    # índice do nó inferior-esquerdo (BL)
            conn[ecount - 1, :] = [n1, n1 + 1, n1 + n_side + 2, n1 + n_side + 1]  # BL BR TR TL
    #----------------------------
    dofs = 3 * n_nos                                           # 3 graus de liberdade por nó (w, theta_x, theta_y)
    K = sp_sparse.lil_matrix((dofs, dofs), dtype=float)        
    F = np.zeros((dofs,), dtype=float)                         # inicializa vetor de forças globais
    #----------------------------
    for e in range(n_quadrados):
        nodes = conn[e, :]                                     # lista de 4 nós do elemento (ainda 1-based)
        # obter coordenadas do elemento para x e y dos 4 nós (presume-se ordem BL BR TR TL)
        Xe = coords[nodes - 1, 0]                              # x coords dos nós do elemento
        Ye = coords[nodes - 1, 1]                              # y coords dos nós do elemento
        x1v = Xe[0]                                            # x1 do elemento (BL.x)
        x2v = Xe[1]                                            # x2 do elemento (BR.x) 
        y1v = Ye[0]                                            # y1 do elemento (BL.y)
        y2v = Ye[3]                                            # y2 do elemento (TL.y)

        I_num_mat = np.array(sp.N(I_xy.subs({x1: x1v, x2: x2v, y1: y1v, y2: y2v})), dtype=float)
        C_num = np.zeros((12, 12), dtype=float)               # inicializa C_num
        C_sub = C_sym.subs({x: x1v, y: y1v})                  # valor simbólico em (x1,y1)  --> BL
        C_num[0:3, :] = np.array(C_sub[0:3, :].evalf(), dtype=float)
        C_sub = C_sym.subs({x: x2v, y: y1v})                  # valor simbólico em (x2,y1)  --> BR
        C_num[3:6, :] = np.array(C_sub[3:6, :].evalf(), dtype=float)
        C_sub = C_sym.subs({x: x2v, y: y2v})                  # valor simbólico em (x2,y2)  --> TR
        C_num[6:9, :] = np.array(C_sub[6:9, :].evalf(), dtype=float)
        C_sub = C_sym.subs({x: x1v, y: y2v})                  # valor simbólico em (x1,y2)  --> TL
        C_num[9:12, :] = np.array(C_sub[9:12, :].evalf(), dtype=float)

        try:
            cond_C = np.linalg.cond(C_num)                     # número de condição
        except np.linalg.LinAlgError:
            cond_C = np.inf
        if cond_C > 1e12:
            warnings.warn(f'C_num mal condicionado no elemento {e+1} (cond={cond_C:.3e}).')
        try:
            Cinv = np.linalg.inv(C_num)                       # inversa de C_num
        except np.linalg.LinAlgError:
            warnings.warn(f'Erro a inverter C_num no elemento {e+1}; será usado pseudo-inversa.')
            Cinv = np.linalg.pinv(C_num)                      #matriz inversa generalizada

        Ke = Cinv.T.dot(I_num_mat).dot(Cinv)                 # montagem de Ke 
        Ae = (x2v - x1v) * (y2v - y1v)                       # área do elemento
        Fe = np.zeros((12,), dtype=float)                    # vetor de cargas elementar
        for ni in range(4):
            Fe[ni*3 + 0] = p0 * Ae / 4.0                     # aplicar carga p0 distribuída igualmente em w-dof

        # espalhar Ke e Fe para K e F globais (localizar os graus de liberdade do elemento)
        dofs_e = np.zeros((12,), dtype=int)                  # índices globais dos graus de liberade do elemento
        for iN in range(4):
            
            node_idx = nodes[iN] - 1                         # converter 1-based -> 0-based
            dofs_e[iN*3 + 0] = node_idx * 3 + 0              # w DOF
            dofs_e[iN*3 + 1] = node_idx * 3 + 1              # theta_x DOF
            dofs_e[iN*3 + 2] = node_idx * 3 + 2              # theta_y DOF

        for r in range(12):
            for c_idx in range(12):
                K[dofs_e[r], dofs_e[c_idx]] += Ke[r, c_idx]  
        for r in range(12):
            F[dofs_e[r]] += Fe[r]

    def apply_borda_fun(bord, coords_local, n_nos_local, dofs_local, tol_local, fixed_in):
        fixed_out = fixed_in.copy()                          # copia do vetor de fixos
        typ = bord['type']                                   # tipo de borda
        coord = bord['coord']                                # coordenada associada ('x' ou 'y')
        val = bord['value']                                  # valor do borda (0 ou a ou b)
        for ni in range(n_nos_local):                        # corre os nós
            x_i = coords_local[ni, 0]                        # coordenada x do nó
            y_i = coords_local[ni, 1]                        # coordenada y do nó
            on_borda = False                                 # flag se nó está na borda
            if coord == 'x':                                 # se bordo vertical (x fixo)
                if abs(x_i - val) < tol_local:
                    on_borda = True
            else:                                            # bordo horizontal (y fixo)
                if abs(y_i - val) < tol_local:
                    on_borda = True
            if on_borda:
                idx_base = ni * 3                            # base dos graus de liberdade do nó no vetor global
                if typ == 'C':                               # encastrado: w, theta_x, theta_y fixos
                    fixed_out[idx_base + 0] = True           # fixar w
                    fixed_out[idx_base + 1] = True           # fixar theta_x
                    fixed_out[idx_base + 2] = True           # fixar theta_y
                elif typ == 'S':                             # simplesmente apoiado
                    if coord == 'y':                         # borda horizontal -> w=0 e theta_x=0
                        fixed_out[idx_base + 0] = True       # w
                        fixed_out[idx_base + 1] = True       # theta_x
                    else:                                    # borda vertical -> w=0 e theta_y=0
                        fixed_out[idx_base + 0] = True       # w
                        fixed_out[idx_base + 2] = True       # theta_y
                elif typ == 'F':                             # livre -> nada a fixar
                    pass
                else:
                    raise ValueError(f"Tipo de borda desconhecido: {typ} (esperado C, S ou F)")
        return fixed_out

    borda_map = {
        'B1': {'type': B1, 'coord': 'x', 'value': 0.0},    # B1: x=0
        'B2': {'type': B2, 'coord': 'x', 'value': a},      # B2: x=a
        'B3': {'type': B3, 'coord': 'y', 'value': 0.0},    # B3: y=0
        'B4': {'type': B4, 'coord': 'y', 'value': b}       # B4: y=b
    }

    tol = 1e-9                                            # tolerância para identificar nós na borda
    fixed = np.zeros((dofs,), dtype=bool)                 
    # aplicar cada borda ao vetor fixed
    fixed = apply_borda_fun(borda_map['B1'], coords, n_nos, dofs, tol, fixed)
    fixed = apply_borda_fun(borda_map['B2'], coords, n_nos, dofs, tol, fixed)
    fixed = apply_borda_fun(borda_map['B3'], coords, n_nos, dofs, tol, fixed)
    fixed = apply_borda_fun(borda_map['B4'], coords, n_nos, dofs, tol, fixed)
    free = np.where(~fixed)[0]                            # índices dos DOFs livres

    K_csr = K.tocsr()                                     # converte para CSR
    Kff = K_csr[free, :][:, free]                         # submatriz livre-livre (esparsa)
    Ff = F[free]                                          # vetor de forças livres

    try:
        U = np.zeros((dofs,), dtype=float)                # inicializar vetor de deslocamentos
        Uf = spla.spsolve(Kff.tocsc(), Ff)                # resolve os graus de liberdade livres
        U[free] = Uf                                      # monta o vetor completo U
    except Exception as ex:
        warnings.warn(f'Falha ao resolver sistema linear esparso: {ex}; tentando solução densa.')
        Kff_dense = Kff.toarray()
        Uf = np.linalg.solve(Kff_dense, Ff)
        U[free] = Uf

    # calcular deslocamentos nodais Wnod (w DOFs, cada 3º valor começando em 0)
    Wnod = U[0::3]                                        # extrair w do vetor U
    w_max = np.max(np.abs(Wnod))                          # valor máximo absoluto da deflexão nodal

    # calcular erro relativo face ao passo anterior (em %)
    if not np.isnan(w_prev):
        erro_rel = abs(w_max - w_prev) / w_prev * 100.0
    else:
        erro_rel = np.inf

    print(f'Iteração: malha {n_side}x{n_side} -> max|w| = {w_max:.6e} m  (erro rel = {erro_rel:.3f}%)')

    if (not np.isnan(w_prev)) and (erro_rel <= tol_rel_pct):
        converged = True
        print(f'\nConvergência atingida! erro relativo = {erro_rel:.3f}%\n')
        # guardar solução convergida para pós-processamento
        U_converged = U.copy()
        coords_converged = coords.copy()
        conn_converged = conn.copy()
        n_nos_converged = n_nos
        K_converged = K_csr.copy()
        F_converged = F.copy()
        Reac_converged = (K_csr.dot(U) - F)               # reações nodais
        break


    w_prev = w_max
# ----------------------------
U = U_converged                                            # vetor de deslocamentos convergido
coords = coords_converged                                  # coordenadas convergidas
conn = conn_converged                                      # conectividade convergida
n_nos = n_nos_converged                                    # número de nós convergido
dofs = 3 * n_nos                                           # graus de liberdade totais
K_full = K_converged                                       # matriz K convergida 
F_full = F_converged                                       # vetor de forças convergido
Reac = Reac_converged                                      # vetor de reações (K*U - F)
Wnod = U[0::3]                                             # deflexões nodais w

Reac_y = Reac[0::3]      # forças verticais (reacções de w)
Reac_tx = Reac[1::3]     # momentos em torno de x (reacções de theta_x)
Reac_ty = Reac[2::3]     # momentos em torno de y (reacções de theta_y)

# Identificar nós encastrados (x=0 e x=a)
tol = 1e-9
n_enc = np.where(
    (np.abs(coords[:,0] - 0.0) < tol) | (np.abs(coords[:,0] - a) < tol)
)[0]

# Extrair apenas as reacções nos nós encastrados
Ry_enc = Reac_y[n_enc]
Mx_enc = Reac_tx[n_enc]
My_enc = Reac_ty[n_enc]

print('\n--- REACÇÕES NOS ENCASTRAMENTOS ---')
print(f'Nº de nós encastrados: {len(n_enc)}')
print(f'Força vertical Ry: max = {np.max(Ry_enc):.3e} N | min = {np.min(Ry_enc):.3e} N')
print(f'Momento Mx: max = {np.max(Mx_enc):.3e} N·m | min = {np.min(Mx_enc):.3e} N·m')
print(f'Momento My: max = {np.max(My_enc):.3e} N·m | min = {np.min(My_enc):.3e} N·m')

# Usar a relção constitutiva para calcular tensões nodais
Qmat = Q_mat
sigma_x_node = np.zeros((n_nos,), dtype=float)              # sigma_x por nó
sigma_y_node = np.zeros((n_nos,), dtype=float)              # sigma_y por nó
tau_xy_node = np.zeros((n_nos,), dtype=float)               # tau_xy por nó
cont = np.zeros((n_nos,), dtype=int)                        # contador de contribuições por nó
gp = np.array([-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)])

# loop por elementos para calcular tensões nodais (média das contribuições)
for e in range(conn.shape[0]):
    nodes = conn[e, :]                                      # nós do elemento (1-based)
    Xe = coords[nodes - 1, 0]                               # x coords dos 4 nós
    Ye = coords[nodes - 1, 1]                               # y coords dos 4 nós
    # deslocamentos nodais do elemento (12x1)
    Ue = np.zeros((12,), dtype=float)
    for iN in range(4):
        node_idx = nodes[iN] - 1
        Ue[iN*3:(iN*3 + 3)] = U[node_idx*3:(node_idx*3 + 3)]
    # montar C_num novamente (como antes)
    x1v = Xe[0]; x2v = Xe[1]; y1v = Ye[0]; y2v = Ye[3]
    C_num = np.zeros((12, 12), dtype=float)
    C_sub = C_sym.subs({x: x1v, y: y1v}); C_num[0:3, :] = np.array(C_sub[0:3, :].evalf(), dtype=float)
    C_sub = C_sym.subs({x: x2v, y: y1v}); C_num[3:6, :] = np.array(C_sub[3:6, :].evalf(), dtype=float)
    C_sub = C_sym.subs({x: x2v, y: y2v}); C_num[6:9, :] = np.array(C_sub[6:9, :].evalf(), dtype=float)
    C_sub = C_sym.subs({x: x1v, y: y2v}); C_num[9:12, :] = np.array(C_sub[9:12, :].evalf(), dtype=float)
    # inversa para obter coeficientes a_elem
    try:
        Cinv = np.linalg.inv(C_num)
    except np.linalg.LinAlgError:
        Cinv = np.linalg.pinv(C_num)
    a_elem = Cinv.dot(Ue)                                   # coeficientes do polinómio do elemento

    # integração 2x2 Gauss para avaliar curvaturas e tensões e acumular em nós
    for i in range(2):
        xi = gp[i]
        for j in range(2):
            eta = gp[j]
            # mapear xi,eta para coordenadas físicas x_p,y_p (assumindo mapeamento bilinear simplificado)
            x_p = 0.5 * ((1 - xi) * x1v + (1 + xi) * x2v)  # mapeamento linear em x (elementos regulares)
            y_p = 0.5 * ((1 - eta) * y1v + (1 + eta) * y2v)  # mapeamento linear em y
            H_eval = np.array(H.subs({x: x_p, y: y_p}).evalf(), dtype=float)  # avalia H no ponto (x_p,y_p)
            kappa = H_eval.dot(a_elem)                       # curvaturas (3x1)
            z = t / 2.0                                      # posição z para cálculo de tensões (face superior = +t/2)
            sigma_vec = -z * (Qmat.dot(kappa))              # vetor [sigma_x; sigma_y; tau_xy] = -z * Q * kappa
            # encontrar nó mais próximo para agregar contributo (simplificação)
            dists = np.sum((coords - np.array([x_p, y_p]))**2, axis=1)
            n_i = np.argmin(dists)                          # índice do nó mais próximo (0-based)
            sigma_x_node[n_i] += sigma_vec[0]               # acumular sigma_x
            sigma_y_node[n_i] += sigma_vec[1]               # acumular sigma_y
            tau_xy_node[n_i] += sigma_vec[2]                # acumular tau_xy
            cont[n_i] += 1                                  # incrementar contador de contributos

# média das contribuições em cada nó (somente onde cont>0)
mask = cont > 0
sigma_x_node[mask] = sigma_x_node[mask] / cont[mask]
sigma_y_node[mask] = sigma_y_node[mask] / cont[mask]
tau_xy_node[mask] = tau_xy_node[mask] / cont[mask]
# ----------------------------
print(f'--- RESULTADOS FINAIS (malha {n_side}x{n_side}) ---')
print(f'Deflexão máxima (|w|) = {np.max(np.abs(Wnod)):.6e} m')
ReacZ = Reac[0:dofs:3]                                       # reações na componente vertical (cada 3º)
print(f'Reacção vertical: max = {np.max(ReacZ):.3e} N | min = {np.min(ReacZ):.3e} N')
if np.any(mask):
    print('Tensões (z = +t/2) (Pa):')
    print(f' sigma_x: max = {np.max(sigma_x_node[mask]):.3e} , min = {np.min(sigma_x_node[mask]):.3e}')
    print(f' sigma_y: max = {np.max(sigma_y_node[mask]):.3e} , min = {np.min(sigma_y_node[mask]):.3e}')
    print(f' tau_xy : max = {np.max(tau_xy_node[mask]):.3e} , min = {np.min(tau_xy_node[mask]):.3e}')
else:
    print('Nenhuma contribuição de tensão calculada (verifica malha).')

# ----------------------------
Xg_lin = np.linspace(0, a, 200)                               # 200 pontos em x
Yg_lin = np.linspace(0, b, 200)                               # 200 pontos em y
Xg, Yg = np.meshgrid(Xg_lin, Yg_lin)                          # grelha 2D
Wgrid = griddata(points=coords, values=Wnod, xi=(Xg, Yg), method='cubic')  

#Gráfico da deflexão
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Xg, Yg, Wgrid, cmap='jet', edgecolor='none')        # superfície
fig.colorbar(surf, ax=ax, shrink=0.6)                                      # barra de cor
ax.set_title('Deflexão w (m) - convergido')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)'); ax.set_zlabel('w (m)')
ax.view_init(elev=45, azim=30)
plt.tight_layout()

#Grafico sigma_x
Sxgrid = griddata(points=coords, values=sigma_x_node, xi=(Xg, Yg), method='cubic')
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
surf2 = ax2.plot_surface(Xg, Yg, Sxgrid, cmap='jet', edgecolor='none')
fig2.colorbar(surf2, ax=ax2, shrink=0.6)
ax2.set_title(r'$\sigma_x$ (Pa) - face superior'); ax2.set_xlabel('x (m)'); ax2.set_ylabel('y (m)')
ax2.view_init(elev=45, azim=30)
plt.tight_layout()

#Grafico sigma_y
Sygrid = griddata(points=coords, values=sigma_y_node, xi=(Xg, Yg), method='cubic')
fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='3d')
surf3 = ax3.plot_surface(Xg, Yg, Sygrid, cmap='jet', edgecolor='none')
fig3.colorbar(surf3, ax=ax3, shrink=0.6)
ax3.set_title(r'$\sigma_y$ (Pa) - face superior'); ax3.set_xlabel('x (m)'); ax3.set_ylabel('y (m)')
ax3.view_init(elev=45, azim=30)
plt.tight_layout()

#Gráfico tau_xy 
Txygrid = griddata(points=coords, values=tau_xy_node, xi=(Xg, Yg), method='cubic')
fig4 = plt.figure()
ax4 = fig4.add_subplot(111, projection='3d')
surf4 = ax4.plot_surface(Xg, Yg, Txygrid, cmap='jet', edgecolor='none')
fig4.colorbar(surf4, ax=ax4, shrink=0.6)
ax4.set_title(r'$\tau_{xy}$ (Pa) - face superior'); ax4.set_xlabel('x (m)'); ax4.set_ylabel('y (m)')
ax4.view_init(elev=45, azim=30)
plt.tight_layout()

#Gráfico Reacções verticais
Rgrid = griddata(points=coords, values=ReacZ, xi=(Xg, Yg), method='cubic')
fig5 = plt.figure()
ax5 = fig5.add_subplot(111, projection='3d')
surf5 = ax5.plot_surface(Xg, Yg, Rgrid, cmap='jet', edgecolor='none')
fig5.colorbar(surf5, ax=ax5, shrink=0.6)
ax5.set_title('Reacções verticais (N)'); ax5.set_xlabel('x (m)'); ax5.set_ylabel('y (m)')
ax5.view_init(elev=45, azim=30)
plt.tight_layout()

#Gráficos dos momentos de reacção
Mxgrid = griddata(points=coords, values=Reac_tx, xi=(Xg, Yg), method='cubic')
Mygrid = griddata(points=coords, values=Reac_ty, xi=(Xg, Yg), method='cubic')

fig6 = plt.figure()
ax6 = fig6.add_subplot(111, projection='3d')
surf6 = ax6.plot_surface(Xg, Yg, Mxgrid, cmap='jet', edgecolor='none')
fig6.colorbar(surf6, ax=ax6, shrink=0.6)
ax6.set_title('Momento de reacção Mx (N·m)')
ax6.set_xlabel('x (m)'); ax6.set_ylabel('y (m)')
ax6.view_init(elev=45, azim=30)
plt.tight_layout()

fig7 = plt.figure()
ax7 = fig7.add_subplot(111, projection='3d')
surf7 = ax7.plot_surface(Xg, Yg, Mygrid, cmap='jet', edgecolor='none')
fig7.colorbar(surf7, ax=ax7, shrink=0.6)
ax7.set_title('Momento de reacção My (N·m)')
ax7.set_xlabel('x (m)'); ax7.set_ylabel('y (m)')
ax7.view_init(elev=45, azim=30)
plt.tight_layout()

plt.show()
