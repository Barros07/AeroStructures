% ==================================
% MÉTODO C/H SIMBÓLICO + CONVERGÊNCIA
%===================================
%Importar bibliotecas
%Definir as porpriedades do material
%Definir a matriz D para o caso ortotrópico
%Definir o tipo de borda em cada lado
%Definição simbolica das variaveis
    %Definir a matriz H simbolica e inicializaçao em 0
    %Definir o polinomio w_e (a1...a12)
    %Definir as rotações theta_x e theta_y (derivadas)
    %construir a matriz simbolica C e inicializaçao em 0
    %Preencher a matriz C simbolica com ([1,4,7,10] -> w; [2,5,8,11] -> theta_x; [3,6,9,12] -> theta_y)
    %Integrar simbolicamente I_xy = ∫∫ H' * D * H dx dy
%Iniciar o metodo de convergencia
    %Definir o erro relativo de convergencia
    %Iniciar a malha com 1 quadrado por lado
    %iniciar o loop de convergencia
        %Incrementar o numero de quadrados por lado
        %Gerar a malha de nos (cada no tem 3 graus de liberdade) e a conectividade
        %Montar a matriz K e o vetor F e inicializa-los em 0
        % Loop pelos elementos para calcular Ke e Fe e espalhar para K, F
        % Obter coordenadas do elemento
        % Avaliar integral simbolica I_xy para o elemento
        % Construir c_num
        % Substituir em cada no e preencher as linhas correspondentes
        % Verificar o condicionamento de C_num e calcular a inversa
        % Calcular a matriz elementar Ke e o vetor de cargas elementar Fe    
        %Aplicar as condições de fronteira
        %Definir a funçao para aplicar as condiçoes de fronteira (seja que tipo de fronteira for)
        %Resolver o sistema de equações lineares 
        %Montar o vetor dos deslocamentos nodais e inicializa em 0 
        %Calcular os deslocamentos nodais w
        %Calcular o erro relativo face ao passo anterior (verificar a convergencia)
%Usar a solução armnazenada
%Calcular as tensões nodais
%Resultados numéricos e impressões finais
%Represenentação gráfica dos resultados

clear; close all; clc;
warning('off','all')

fprintf('=== CASO 6 (método C/H simbólico, convergente) ===\n');

% ----------------------------
a = 0.750;          % [m]
b = 0.600;          % [m]
t = 0.002;          % [m]
E1 = 130e9;         % [Pa]
E2 = 10e9;          % [Pa]
nu12 = 0.26;
G12 = 5e9;          % [Pa]
p0 = -175.0;        % [N/m^2] (negativa = para baixo)
% ----------------------------
nu21 = (E2 / E1) * nu12;
Q11 = E1 / (1 - nu12 * nu21);
Q22 = E2 / (1 - nu12 * nu21);
Q12 = nu21 * Q11;
Q66 = G12;
Q_mat = [Q11, Q12, 0.0;
         Q12, Q22, 0.0;
         0.0,  0.0, Q66];
D_b = (t^3) / 12.0 * Q_mat;
fprintf('\nMatriz D_b (flexão) [Pa·m³]:\n');
disp(D_b);
% ----------------------------
B1 = 'C';   % borda x=0. encastrada
B2 = 'C';   % borda x=a. encastrada
B3 = 'F';   % borda y=0. livre
B4 = 'F';   % borda y=b. livre
% ----------------------------
% Símbolos
syms x y x1 x2 y1 y2 real
D_sym = sym(D_b);

% ----------------------------
H = sym(zeros(3,12));
H(1,:) = [0, 0, 0, -2, 0, 0, -6*x, -2*y, 0, 0, -6*x*y, 0];
H(2,:) = [0, 0, 0, 0, 0, -2, 0, 0, -2*x, -6*y, 0, -6*x*y];
H(3,:) = [0, 0, 0, 0, -2, 0, 0, -4*x, -4*y, 0, -12*x, -12*y]; % notei -6*x*2 -> -12*x
%----------------------------
a_syms = sym('a1:12'); % a1..a12
% construir w_e
w_e = a_syms(1) + a_syms(2)*x + a_syms(3)*y + a_syms(4)*x^2 + a_syms(5)*x*y + a_syms(6)*y^2 + ...
      a_syms(7)*x^3 + a_syms(8)*x^2*y + a_syms(9)*x*y^2 + a_syms(10)*y^3 + ...
      a_syms(11)*x^3*y + a_syms(12)*x*y^3;
%----------------------------
theta_x = diff(w_e, x);    % theta_x = ∂w/∂x
theta_y = diff(w_e, y);    % theta_y = ∂w/∂y
%----------------------------
C_sym = sym(zeros(12,12));
for row = 1:12
    for col = 1:12
        if ismember(row, [1,4,7,10])
            C_sym(row,col) = diff(w_e, a_syms(col));
        elseif ismember(row, [2,5,8,11])
            C_sym(row,col) = diff(theta_x, a_syms(col));
        else
            C_sym(row,col) = diff(theta_y, a_syms(col));
        end
    end
end

%----------------------------
fprintf('A integrar simbolicamente H''*D*H (apenas uma vez) ...\n');
Integrand = H.' * D_sym * H;
I_xy = int(int(Integrand, x, x1, x2), y, y1, y2);
I_xy = simplify(I_xy);

% ----------------------------
tol_rel_pct = 0.5;    % tolerância relativa em percentagem (0.5%)
w_prev = NaN;         % valor anterior w_max (para comparar)
converged = false;
n_side = 1;           % inicia com 1 quadrado por lado (malha 1x1)
fprintf('\nIniciando sequência de refinamento automática (critério: erro <= %.2f%%)\n\n', tol_rel_pct);

while ~converged
    n_side = n_side + 1;                     % incrementa número de quadrados por lado
    n_quadrados = n_side * n_side;
    nnx = n_side + 1;
    nny = n_side + 1;
    n_nos = nnx * nny;
    dx = a / n_side;
    dy = b / n_side;

    coords = zeros(n_nos,2);
    idx = 1;
    for j = 0:n_side
        for i = 0:n_side
            coords(idx,:) = [i*dx, j*dy];
            idx = idx + 1;
        end
    end

    conn = zeros(n_quadrados,4);
    ecount = 0;
    for j = 1:n_side
        for i = 1:n_side
            ecount = ecount + 1;
            n1 = (j-1)*(n_side+1) + i;
            conn(ecount,:) = [n1, n1+1, n1+n_side+2, n1+n_side+1]; % BL BR TR TL (1-based)
        end
    end
    %----------------------------
    dofs = 3 * n_nos;
    K = sparse(dofs, dofs);    % global stiffness (sparse)
    F = zeros(dofs,1);         % global force

    %----------------------------
    for e = 1:n_quadrados
        nodes = conn(e,:);
        Xe = coords(nodes,1);
        Ye = coords(nodes,2);
        x1v = Xe(1);
        x2v = Xe(2);
        y1v = Ye(1);
        y2v = Ye(4);

        % avaliar integral simbolica para o elemento (substituir limites)
        I_num_mat_sym = subs(I_xy, [x1,x2,y1,y2], [x1v, x2v, y1v, y2v]);
        I_num_mat = double(I_num_mat_sym);

        % construir C_num (12x12) usando C_sym avaliado nos 4 nós
        C_num = zeros(12,12);
        C_sub = subs(C_sym, [x,y], [x1v, y1v]);
        C_num(1:3,:) = double(C_sub(1:3,:));
        C_sub = subs(C_sym, [x,y], [x2v, y1v]);
        C_num(4:6,:) = double(C_sub(4:6,:));
        C_sub = subs(C_sym, [x,y], [x2v, y2v]);
        C_num(7:9,:) = double(C_sub(7:9,:));
        C_sub = subs(C_sym, [x,y], [x1v, y2v]);
        C_num(10:12,:) = double(C_sub(10:12,:));

        % verificar condicionamento
        try
            cond_C = cond(C_num);
        catch
            cond_C = Inf;
        end
        if cond_C > 1e12
            warning('C_num mal condicionado no elemento %d (cond=%.3e).', e, cond_C);
        end
        % inversa ou pseudo-inversa
        if rcond(C_num) < 1e-15
            Cinv = pinv(C_num);
            warning('Usada pseudo-inversa para C_num no elemento %d', e);
        else
            Cinv = inv(C_num);
        end

        Ke = Cinv.' * I_num_mat * Cinv;   % matriz elementar
        Ae = (x2v - x1v) * (y2v - y1v);
        Fe = zeros(12,1);
        for ni = 1:4
            Fe((ni-1)*3 + 1) = p0 * Ae / 4.0;  % aplicar carga p0 distribuída igualmente em w-dof
        end

        % espalhar Ke e Fe para K e F globais
        dofs_e = zeros(12,1);
        for iN = 1:4
            node_idx = nodes(iN); % já 1-based
            dofs_e((iN-1)*3 + 1) = (node_idx-1)*3 + 1; % w DOF
            dofs_e((iN-1)*3 + 2) = (node_idx-1)*3 + 2; % theta_x DOF
            dofs_e((iN-1)*3 + 3) = (node_idx-1)*3 + 3; % theta_y DOF
        end

        % assemble (note: MATLAB sparse supports incremental assignment)
        for r = 1:12
            for c_idx = 1:12
                K(dofs_e(r), dofs_e(c_idx)) = K(dofs_e(r), dofs_e(c_idx)) + Ke(r,c_idx);
            end
        end
        for r = 1:12
            F(dofs_e(r)) = F(dofs_e(r)) + Fe(r);
        end
    end

    % funçao para aplicar tipo de borda (retorna vetor lógico fixed)
    function fixed_out = apply_borda_fun_matlab(bord, coords_local, n_nos_local, tol_local, fixed_in)
        fixed_out = fixed_in;
        typ = bord.type;
        coord = bord.coord;
        val = bord.value;
        for ni = 1:n_nos_local
            x_i = coords_local(ni,1);
            y_i = coords_local(ni,2);
            on_borda = false;
            if strcmp(coord,'x')
                if abs(x_i - val) < tol_local
                    on_borda = true;
                end
            else
                if abs(y_i - val) < tol_local
                    on_borda = true;
                end
            end
            if on_borda
                idx_base = (ni-1)*3;
                if strcmp(typ,'C')
                    fixed_out(idx_base + 1) = true;
                    fixed_out(idx_base + 2) = true;
                    fixed_out(idx_base + 3) = true;
                elseif strcmp(typ,'S')
                    if strcmp(coord,'y')
                        fixed_out(idx_base + 1) = true;
                        fixed_out(idx_base + 2) = true;
                    else
                        fixed_out(idx_base + 1) = true;
                        fixed_out(idx_base + 3) = true;
                    end
                elseif strcmp(typ,'F')
                    % nada a fixar
                else
                    error('Tipo de borda desconhecido: %s (esperado C, S ou F)', typ);
                end
            end
        end
    end

    borda_map.B1.type = B1; borda_map.B1.coord = 'x'; borda_map.B1.value = 0.0;
    borda_map.B2.type = B2; borda_map.B2.coord = 'x'; borda_map.B2.value = a;
    borda_map.B3.type = B3; borda_map.B3.coord = 'y'; borda_map.B3.value = 0.0;
    borda_map.B4.type = B4; borda_map.B4.coord = 'y'; borda_map.B4.value = b;

    tol = 1e-9;
    fixed = false(dofs,1);
    fixed = apply_borda_fun_matlab(borda_map.B1, coords, n_nos, tol, fixed);
    fixed = apply_borda_fun_matlab(borda_map.B2, coords, n_nos, tol, fixed);
    fixed = apply_borda_fun_matlab(borda_map.B3, coords, n_nos, tol, fixed);
    fixed = apply_borda_fun_matlab(borda_map.B4, coords, n_nos, tol, fixed);
    free = find(~fixed);

    % submatriz livre-livre
    Kff = K(free, free);
    Ff = F(free);

    % resolver sistema (sparse preferencialmente)
    U = zeros(dofs,1);
    try
        Uf = Kff \ Ff;
        U(free) = Uf;
    catch ME
        warning('Falha ao resolver sistema esparso: %s. Tentando denso...', ME.message);
        Kff_dense = full(Kff);
        Uf = Kff_dense \ Ff;
        U(free) = Uf;
    end

    % calcular deslocamentos nodais Wnod
    Wnod = U(1:3:end);
    w_max = max(abs(Wnod));

    if ~isnan(w_prev)
        erro_rel = abs(w_max - w_prev) / w_prev * 100.0;
    else
        erro_rel = Inf;
    end

    fprintf('Iteração: malha %dx%d -> max|w| = %.6e m  (erro rel = %.3f%%)\n', n_side, n_side, w_max, erro_rel);

    if (~isnan(w_prev)) && (erro_rel <= tol_rel_pct)
        converged = true;
        fprintf('\nConvergência atingida! erro relativo = %.3f%%\n\n', erro_rel);
        U_converged = U;
        coords_converged = coords;
        conn_converged = conn;
        n_nos_converged = n_nos;
        K_converged = K;
        F_converged = F;
        Reac_converged = K * U - F;   % reações nodais
        break;
    end

    w_prev = w_max;
end

% ----------------------------
U = U_converged;
coords = coords_converged;
conn = conn_converged;
n_nos = n_nos_converged;
dofs = 3 * n_nos;
K_full = K_converged;
F_full = F_converged;
Reac = Reac_converged;
Wnod = U(1:3:end);

% Reacções por tipo
Reac_y = Reac(1:3:end);   % forças verticais (w)
Reac_tx = Reac(2:3:end);  % momentos em torno de x (theta_x)
Reac_ty = Reac(3:3:end);  % momentos em torno de y (theta_y)

% Identificar nós encastrados (x=0 e x=a)
tol = 1e-9;
n_enc = find( (abs(coords(:,1) - 0.0) < tol) | (abs(coords(:,1) - a) < tol) );

Ry_enc = Reac_y(n_enc);
Mx_enc = Reac_tx(n_enc);
My_enc = Reac_ty(n_enc);

fprintf('\n--- REACÇÕES NOS ENCASTRAMENTOS ---\n');
fprintf('Nº de nós encastrados: %d\n', length(n_enc));
fprintf('Força vertical Ry: max = %.3e N | min = %.3e N\n', max(Ry_enc), min(Ry_enc));
fprintf('Momento Mx: max = %.3e N·m | min = %.3e N·m\n', max(Mx_enc), min(Mx_enc));
fprintf('Momento My: max = %.3e N·m | min = %.3e N·m\n', max(My_enc), min(My_enc));

% Usar a relção constitutiva para calcular tensões nodais
Qmat = Q_mat;
sigma_x_node = zeros(n_nos,1);
sigma_y_node = zeros(n_nos,1);
tau_xy_node  = zeros(n_nos,1);
cont = zeros(n_nos,1);
gp = [-1.0/sqrt(3.0), 1.0/sqrt(3.0)];

% loop por elementos para calcular tensões nodais (média das contribuições)
for e = 1:size(conn,1)
    nodes = conn(e,:);
    Xe = coords(nodes,1);
    Ye = coords(nodes,2);
    % deslocamentos nodais do elemento (12x1)
    Ue = zeros(12,1);
    for iN = 1:4
        node_idx = nodes(iN);
        base = (iN-1)*3;
        % extrair do vetor global U (note offsets)
        Ue(base+1:base+3) = U((node_idx-1)*3 + (1:3));
    end
    % montar C_num novamente (como antes)
    x1v = Xe(1); x2v = Xe(2); y1v = Ye(1); y2v = Ye(4);
    C_num = zeros(12,12);
    C_sub = subs(C_sym, [x,y], [x1v, y1v]); C_num(1:3,:) = double(C_sub(1:3,:));
    C_sub = subs(C_sym, [x,y], [x2v, y1v]); C_num(4:6,:) = double(C_sub(4:6,:));
    C_sub = subs(C_sym, [x,y], [x2v, y2v]); C_num(7:9,:) = double(C_sub(7:9,:));
    C_sub = subs(C_sym, [x,y], [x1v, y2v]); C_num(10:12,:) = double(C_sub(10:12,:));
    if rcond(C_num) < 1e-15
        Cinv = pinv(C_num);
    else
        Cinv = inv(C_num);
    end
    a_elem = Cinv * Ue;  % coeficientes do polinómio do elemento

    % integração 2x2 Gauss para avaliar curvaturas e tensões e acumular em nós
    for i = 1:2
        xi = gp(i);
        for j = 1:2
            eta = gp(j);
            x_p = 0.5 * ((1 - xi) * x1v + (1 + xi) * x2v);  % mapeamento linear em x
            y_p = 0.5 * ((1 - eta) * y1v + (1 + eta) * y2v);% mapeamento linear em y
            H_eval_sym = subs(H, [x,y], [x_p, y_p]);
            H_eval = double(H_eval_sym);
            kappa = H_eval * a_elem; % curvaturas (3x1)
            z = t / 2.0;
            sigma_vec = -z * (Qmat * kappa); % [sigma_x; sigma_y; tau_xy]
            % encontrar nó mais próximo para agregar contributo (simplificação)
            dists = sum((coords - repmat([x_p, y_p], n_nos, 1)).^2, 2);
            [~, n_i] = min(dists);
            sigma_x_node(n_i) = sigma_x_node(n_i) + sigma_vec(1);
            sigma_y_node(n_i) = sigma_y_node(n_i) + sigma_vec(2);
            tau_xy_node(n_i)  = tau_xy_node(n_i)  + sigma_vec(3);
            cont(n_i) = cont(n_i) + 1;
        end
    end
end

% média das contribuições em cada nó (somente onde cont>0)
mask = cont > 0;
sigma_x_node(mask) = sigma_x_node(mask) ./ cont(mask);
sigma_y_node(mask) = sigma_y_node(mask) ./ cont(mask);
tau_xy_node(mask)  = tau_xy_node(mask)  ./ cont(mask);

% ----------------------------
fprintf('--- RESULTADOS FINAIS (malha %dx%d) ---\n', n_side, n_side);
fprintf('Deflexão máxima (|w|) = %.6e m\n', max(abs(Wnod)));
ReacZ = Reac(1:3:end);
fprintf('Reacção vertical: max = %.3e N | min = %.3e N\n', max(ReacZ), min(ReacZ));
if any(mask)
    fprintf('Tensões (z = +t/2) (Pa):\n');
    fprintf(' sigma_x: max = %.3e , min = %.3e\n', max(sigma_x_node(mask)), min(sigma_x_node(mask)));
    fprintf(' sigma_y: max = %.3e , min = %.3e\n', max(sigma_y_node(mask)), min(sigma_y_node(mask)));
    fprintf(' tau_xy : max = %.3e , min = %.3e\n', max(tau_xy_node(mask)), min(tau_xy_node(mask)));
else
    fprintf('Nenhuma contribuição de tensão calculada (verifica malha).\n');
end

% ----------------------------
Xg_lin = linspace(0, a, 200);
Yg_lin = linspace(0, b, 200);
[Xg, Yg] = meshgrid(Xg_lin, Yg_lin);
Wgrid = griddata(coords(:,1), coords(:,2), Wnod, Xg, Yg, 'cubic');

%Gráfico da deflexão
figure;
surf(Xg, Yg, Wgrid, 'EdgeColor', 'none');
colormap('jet'); colorbar;
title('Deflexão w (m) - convergido');
xlabel('x (m)'); ylabel('y (m)'); zlabel('w (m)');
view(45,30); axis tight;

%Grafico sigma_x
Sxgrid = griddata(coords(:,1), coords(:,2), sigma_x_node, Xg, Yg, 'cubic');
figure;
surf(Xg, Yg, Sxgrid, 'EdgeColor', 'none'); colorbar;
title('\sigma_x (Pa) - face superior'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

%Grafico sigma_y
Sygrid = griddata(coords(:,1), coords(:,2), sigma_y_node, Xg, Yg, 'cubic');
figure;
surf(Xg, Yg, Sygrid, 'EdgeColor', 'none'); colorbar;
title('\sigma_y (Pa) - face superior'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

%Gráfico tau_xy
Txygrid = griddata(coords(:,1), coords(:,2), tau_xy_node, Xg, Yg, 'cubic');
figure;
surf(Xg, Yg, Txygrid, 'EdgeColor', 'none'); colorbar;
title('\tau_{xy} (Pa) - face superior'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

%Gráfico Reacções verticais
Rgrid = griddata(coords(:,1), coords(:,2), ReacZ, Xg, Yg, 'cubic');
figure;
surf(Xg, Yg, Rgrid, 'EdgeColor', 'none'); colorbar;
title('Reacções verticais (N)'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

%Gráficos dos momentos de reacção
Mxgrid = griddata(coords(:,1), coords(:,2), Reac_tx, Xg, Yg, 'cubic');
Mygrid = griddata(coords(:,1), coords(:,2), Reac_ty, Xg, Yg, 'cubic');

figure;
surf(Xg, Yg, Mxgrid, 'EdgeColor', 'none'); colorbar;
title('Momento de reacção Mx (N·m)'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

figure;
surf(Xg, Yg, Mygrid, 'EdgeColor', 'none'); colorbar;
title('Momento de reacção My (N·m)'); xlabel('x (m)'); ylabel('y (m)');
view(45,30); axis tight;

% show done
fprintf('Processo terminado.\n');
