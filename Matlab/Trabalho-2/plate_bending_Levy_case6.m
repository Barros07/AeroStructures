% plate_bending_Levy_case6.m
% Versão MATLAB Online compatível, com soluções simbólicas e exportação para ficheiro

%% 0) Limpar e definir pacotes
clear; close all; clc;

%% 1) Símbolos e solução simbólica por modo
syms Am Bm Cm Dm x y b a alpha nu D Pm kw ktheta real

% Solução particular (constante em y)
W_part = (Pm / D) / alpha^4;

% Solução homogénea
W_h = Am*cosh(alpha*y) + Bm*sinh(alpha*y) + Cm*y.*cosh(alpha*y) + Dm*y.*sinh(alpha*y);

% Solução total
W_total = W_h + W_part;

% Deflexão completa: w(x,y) = W_total(y) * sin(alpha*x)
w_xy = W_total * sin(alpha*x);

%% 1.3) Derivadas necessárias
w_y   = diff(w_xy, y);
w_yy  = diff(w_y, y);
w_yyy = diff(w_yy, y);
w_xx  = diff(w_xy, x, 2);

%% 1.4) Avaliar expressões em y = 0 e y = b
w0     = simplify(subs(w_xy, y, 0));
w_y0   = simplify(subs(w_y, y, 0));
w_yy0  = simplify(subs(w_yy, y, 0));
w_yyy0 = simplify(subs(w_yyy, y, 0));

w_b     = simplify(subs(w_xy, y, b));
w_yb   = simplify(subs(w_y, y, b));
w_yyb  = simplify(subs(w_yy, y, b));
w_yyyb = simplify(subs(w_yyy, y, b));

%% 1.5) Momentos e esforço cortante
M_y0 = D*(w_yy0 - nu*subs(w_xx, y, 0));
M_yb = D*(w_yyb - nu*subs(w_xx, y, b));

Q_y0 = D*(w_yyy0 - (2 - nu)*alpha^2 * w_y0);
Q_yb = D*(w_yyyb - (2 - nu)*alpha^2 * w_yb);

%% 1.6) Condições de fronteira
eq1_x = simplify(Q_y0 + kw * w0);
eq2_x = simplify(M_y0);
eq3_x = simplify(Q_yb - kw * w_b);
eq4_x = simplify(M_yb - ktheta * w_yb);

sin_factor = sin(alpha*x);
eq1 = simplify(eq1_x / sin_factor);
eq2 = simplify(eq2_x / sin_factor);
eq3 = simplify(eq3_x / sin_factor);
eq4 = simplify(eq4_x / sin_factor);

%% 1.7) Montar matriz simbólica e vetor RHS
unknowns = [Am, Bm, Cm, Dm];
eqs = {eq1, eq2, eq3, eq4};
M_sym = sym(zeros(4,4));
RHS_sym = sym(zeros(4,1));

for i = 1:4
    eq = expand(eqs{i});
    for j = 1:4
        % Tentar extrair coeficiente linear de cada incógnita
        C_all = coeffs(eq, unknowns(j));
        if numel(C_all) > 1
            M_sym(i,j) = C_all(2);
        else
            if has(eq, unknowns(j))
                M_sym(i,j) = C_all;
            else
                M_sym(i,j) = sym(0);
            end
        end
    end
    const_term = subs(eq, {Am,Bm,Cm,Dm}, {0,0,0,0});
    RHS_sym(i) = -const_term;
end

%% 1.8) Solução simbólica geral
sol_struct = solve(M_sym * [Am;Bm;Cm;Dm] == RHS_sym, [Am,Bm,Cm,Dm], 'IgnoreAnalyticConstraints', true);

A_expr = simplifyFraction(simplify(sol_struct.Am, 'Steps', 50));
B_expr = simplifyFraction(simplify(sol_struct.Bm, 'Steps', 50));
C_expr = simplifyFraction(simplify(sol_struct.Cm, 'Steps', 50));
D_expr = simplifyFraction(simplify(sol_struct.Dm, 'Steps', 50));

% Organizar por parâmetro alpha (melhora legibilidade)
A_expr = collect(A_expr, alpha);
B_expr = collect(B_expr, alpha);
C_expr = collect(C_expr, alpha);
D_expr = collect(D_expr, alpha);

disp('--- Soluções simbólicas gerais (simplificadas) ---');
fprintf('\nA_m = %s\n', char(A_expr));
fprintf('B_m = %s\n', char(B_expr));
fprintf('C_m = %s\n', char(C_expr));
fprintf('D_m = %s\n\n', char(D_expr));

%% 2) Valores do Caso 6 (numérico)
a_val = 0.75; b_val = 0.6; t_val = 0.001;
E_val = 210e9; nu_val = 0.30; kw_val = 125e3; ktheta_val = 1060; p0_val = 175;
D_val = E_val * t_val^3 / (12*(1 - nu_val^2));

compute_Pm_alpha = @(m,p0,a) deal(4.0 * p0 / (a * (m*pi/a)), m*pi / a);

m_list = [1,3,5,7,9];
ABCD = struct();

for idx = 1:length(m_list)
    m_val = m_list(idx);
    [Pm_val, alpha_val] = compute_Pm_alpha(m_val, p0_val, a_val);

    subs_map = [alpha, a, b, D, nu, kw, ktheta, Pm];
    subs_vals = [alpha_val, a_val, b_val, D_val, nu_val, kw_val, ktheta_val, Pm_val];

    M_num_sym = subs(M_sym, subs_map, subs_vals);
    RHS_num_sym = subs(RHS_sym, subs_map, subs_vals);

    M_num = double(M_num_sym);
    RHS_num = double(RHS_num_sym);

    if rcond(M_num) > eps
        sol = M_num \ RHS_num;
    else
        sol = pinv(M_num) * RHS_num;
    end

    ABCD(m_val).A = sol(1);
    ABCD(m_val).B = sol(2);
    ABCD(m_val).C = sol(3);
    ABCD(m_val).D = sol(4);
    ABCD(m_val).Pm = Pm_val;
    ABCD(m_val).alpha = alpha_val;
end

%% 3) Soma modal e plot
nx = 100; ny = 100;
x_vals = linspace(0, a_val, nx);
y_vals = linspace(0, b_val, ny);
[X, Y] = meshgrid(x_vals, y_vals);
W = zeros(size(X));

for m_val = m_list
    s = ABCD(m_val);
    alpha_m = s.alpha; Pm_m = s.Pm; A = s.A; B = s.B; C = s.C; Dc = s.D;
    cosh_term = cosh(alpha_m .* Y);
    sinh_term = sinh(alpha_m .* Y);
    Wm = (A .* cosh_term + B .* sinh_term + C .* Y .* cosh_term + Dc .* Y .* sinh_term + Pm_m ./ (D_val * alpha_m^4));
    W = W + Wm .* sin(alpha_m .* X);
end

W = -W;
figure('Name','Deflexão w(x,y) - Caso 6 (Lévy) - Corrigido','NumberTitle','off')
surf(X,Y,W,'EdgeColor','none')
colormap turbo; colorbar
xlabel('x (m)'); ylabel('y (m)'); zlabel('w (m)')
title('Deflexão w(x,y) - Caso 6 (Lévy) - Corrigido')
view(3)

[~, ix] = min(abs(x_vals - a_val/2));
[~, iy] = min(abs(y_vals - b_val/2));
w_centro = W(iy, ix);
fprintf('w_centro ≈ %.3e m (negativo = para baixo)\n', w_centro);

%% 4) Exportar resultados para ficheiro de texto
fname = 'resultados_Levy_Caso6.txt';
fid = fopen(fname, 'w');

fprintf(fid, 'Resultados - Plate Bending (Lévy) - Caso 6\n');
fprintf(fid, '------------------------------------------\n');
fprintf(fid, 'a = %.3f m, b = %.3f m, t = %.4f m\n', a_val, b_val, t_val);
fprintf(fid, 'E = %.2e Pa, nu = %.2f\n', E_val, nu_val);
fprintf(fid, 'kw = %.2e N/m^3, ktheta = %.2f N·m/rad, p0 = %.2f N/m^2\n', kw_val, ktheta_val, p0_val);
fprintf(fid, 'D = %.3e N·m\n\n', D_val);

fprintf(fid, '--- Soluções simbólicas gerais ---\n');
fprintf(fid, 'A_m = %s\n', char(A_expr));
fprintf(fid, 'B_m = %s\n', char(B_expr));
fprintf(fid, 'C_m = %s\n', char(C_expr));
fprintf(fid, 'D_m = %s\n\n', char(D_expr));

fprintf(fid, 'Constantes por modo:\n');
fprintf(fid, ' m |       A [m]        |       B [m]        |       C [m]        |       D [m]\n');
fprintf(fid, '---|--------------------|--------------------|--------------------|--------------------\n');
for m_val = m_list
    s = ABCD(m_val);
    fprintf(fid, '%2d | %+.5e | %+.5e | %+.5e | %+.5e\n', m_val, s.A, s.B, s.C, s.D);
end

fprintf(fid, '\nDeflexão no centro:\n');
fprintf(fid, 'w_centro = %.5e m (negativo = para baixo)\n', w_centro);
fclose(fid);

disp(['Resultados exportados para "', fname, '" com sucesso.']);
