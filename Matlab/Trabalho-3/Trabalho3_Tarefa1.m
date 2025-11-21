% ==================================
% MÉTODO DE RITZ — PLACA ORTOTRÓPICA
% ==================================

%Introduçao de bibliotecas
clear; clc; close all;

%Introdução dos parâmetros geométricos 
%Calculo das rigidezes de flexão
%Pedir o numero de modos
%Definição da função X e Y
%Definiçao da funçao do deslocamento w
%Definição das matrizes R e B
%Resolver o problema generalizado de autovalores
%Pedir ao utilzador o modo a visualizar
%Reconstruir a funçao w(x,y)
%Reprensentação simbólica da frequência ω
%Cálculo do deslocamento máximo |w|max
%Representação gráfica do modo de vibração


% =========================
% 1) PARÂMETROS GEOMÉTRICOS
% =========================
a = 0.750;     % [m]
b = 0.600;     % [m]
t = 0.002;     % [m]
rho = 1600;    % [kg/m³]
E1 = 130e9;    % [Pa]
E2 = 10e9;     % [Pa]
nu12 = 0.26; 
G12 = 5e9;     % [Pa]

nu21 = nu12 * E2 / E1;              

Q11 = E1 / (1 - nu12 * nu21);
Q22 = E2 / (1 - nu12 * nu21);
Q12 = nu12 * E2 / (1 - nu12 * nu21);
Q66 = G12;

D11 = Q11 * t^3 / 12;
D22 = Q22 * t^3 / 12;
D12 = Q12 * t^3 / 12;
D66 = Q66 * t^3 / 12;


% ========================
% 2) PEDIR NÚMERO DE MODOS
% ========================
n_modos = input('Introduz o número de modos de vibração: ');
if isempty(n_modos) || n_modos < 1
    n_modos = 2;
end

M = n_modos;
N = n_modos;
n = M * N;   % número total de graus de liberdade


% ============================
% 3) Definição da função X e Y
% ============================
syms x y
X = sym(zeros(1, M));
Y = sym(zeros(1, N));

for i = 1:M
    X(i) = (x / a)^(i + 1) - 2 * (x / a)^(i + 2) + (x / a)^(i + 3);
end

for j = 1:N
    Y(j) = (y / b)^j - (y / b)^(j + 1);
end


% ======================================
% 4) FUNÇÃO DE RITZ (w = Σ c_ij X_i Y_j)
% ======================================
c = sym('c', [n, 1]);
w = sym(0);
k = 1;

for i = 1:M
    for j = 1:N
        w = w + c(k) * X(i) * Y(j);
        k = k + 1;
    end
end


% ============================
% 5) MATRIZES [R] e [B] (Ritz)
% ============================
R = sym(zeros(n));
B = sym(zeros(n));

for p = 1:n
    wp = diff(w, c(p));
    wp_xx = diff(wp, x, 2);
    wp_yy = diff(wp, y, 2);
    wp_xy = diff(diff(wp, x), y);

    for q = 1:n
        wq = diff(w, c(q));
        wq_xx = diff(wq, x, 2);
        wq_yy = diff(wq, y, 2);
        wq_xy = diff(diff(wq, x), y);

        expr_R = D11 * wp_xx * wq_xx + 2 * D12 * wp_xx * wq_yy + ...
                 D22 * wp_yy * wq_yy + 4 * D66 * wp_xy * wq_xy;

        expr_B = rho * t * wp * wq;

        R(p, q) = int(int(expr_R, x, 0, a), y, 0, b);
        B(p, q) = int(int(expr_B, x, 0, a), y, 0, b);
    end
end

R = double(vpa(R, 6));
B = double(vpa(B, 6));


% ===============================================
% 6) RESOLVER PROBLEMA GENERALIZADO [R]c = ω²[B]c
% ===============================================
[vec, val] = eig(R, B);
omega2 = diag(val);
omega2 = real(omega2(omega2 > 1e-8));
omega = sqrt(omega2);
freq = omega / (2 * pi);

fprintf('\n=== FREQUÊNCIAS NATURAIS ===\n');
for k = 1:length(freq)
    fprintf('Modo %d: f = %.3f Hz\n', k, freq(k));
end


% =============================
% 7) ESCOLHER MODO A VISUALIZAR
% =============================
modo = input('\nIntroduz o número do modo a visualizar (ex: 1): ');
if isempty(modo) || modo < 1 || modo > length(freq)
    modo = 1;
end

phi = vec(:, modo);
phi = phi / max(abs(phi));


% =====================
% 8) RECONSTRUIR w(x,y)
% =====================
Nx = 80; Ny = 80;
x_vals = linspace(0, a, Nx);
y_vals = linspace(0, b, Ny);
[Xgrid, Ygrid] = meshgrid(x_vals, y_vals);
W = zeros(Ny, Nx);

k = 1;
for i = 1:M
    for j = 1:N
        X_fun = matlabFunction(X(i), 'Vars', x);
        Y_fun = matlabFunction(Y(j), 'Vars', y);
        W = W + phi(k) * X_fun(Xgrid) .* Y_fun(Ygrid);
        k = k + 1;
    end
end


% ==========================================================
% 9) EXPRESSÃO SIMBÓLICA DA FREQUÊNCIA ω (SEM SUBSTITUIR VALORES)
% ==========================================================

% Define novamente as variáveis simbólicas (sem valores numéricos)
syms x y c1 a b t rho D11 D22 D12 D66 real

% Define as funções de forma genéricas
X1 = (x / a)^2 - 2 * (x / a)^3 + (x / a)^4;
Y1 = (y / b) - (y / b)^2;
w_sym = c1 * X1 * Y1;

% Derivadas de segunda ordem
w_xx = diff(w_sym, x, 2);
w_yy = diff(w_sym, y, 2);
w_xy = diff(diff(w_sym, x), y);

% Energias de deformação e cinética (mantendo símbolos)
U = int(int(D11*w_xx^2 + 2*D12*w_xx*w_yy + D22*w_yy^2 + 4*D66*w_xy^2, y, 0, b), x, 0, a);
T = int(int(rho * t * w_sym^2, y, 0, b), x, 0, a);

% Expressão simbólica da frequência angular
omega_symbolic = simplify(sqrt(U / T));

disp(' ');
disp('=== EXPRESSÃO SIMBÓLICA DA FREQUÊNCIA ω ===');
pretty(omega_symbolic);
disp(' ');
disp('=== EXPRESSÃO EM LATEX ===');
disp(latex(omega_symbolic));


% ==========================================================
% 10) CÁLCULO DO DESLOCAMENTO MÁXIMO |w|max
% ==========================================================
w_max = max(max(abs(W)));
fprintf('\nValor máximo do deslocamento relativo |w|max = %.4f\n', w_max);


% ==========================================================
% 11) VISUALIZAÇÃO GRÁFICA (FUNDO PRETO)
% ==========================================================
figure('Name','Modo de vibração','NumberTitle','off','Color','k');
surf(Xgrid, Ygrid, W, 'EdgeColor', 'none');
colormap turbo
xlabel('x (m)', 'Color', 'w');
ylabel('y (m)', 'Color', 'w');
zlabel('w(x,y)', 'Color', 'w');
title(sprintf('Modo %d — f = %.3f Hz  |w|max = %.4f', modo, freq(modo), w_max), 'Color', 'w');
ax = gca;
ax.Color = 'k';         
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
shading interp;
colorbar;
axis tight;
view(45, 30);
