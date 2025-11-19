clc;
clear all;
% ========================================================== %
% Objetivo: Resolver o sistema de equações lineares pelos métodos e definir a composição da corrente de saída:
%           1) Método de Cramer;
%           2) Método de Eliminação Gaussiana;
%           3) Método de Matriz inversa.
% ========================================================== %

% ======================== Reações ========================= %
% Adotando uma série de reações no processo:
% Reação 1: a*C3H8 + b1*O2 ----> c*C02 + d*H2O
%          [   1       5           3       4  ]
% Reação 2: a*C3H8 + b2*O2 ----> c*C0 + d*H2O
%          [   1      3.5          3      4   ]*(0.324/2) --> Estequiometria
% Reação 3: a*C3H8 + b3*O2 ----> c*C0 + d*H2
%          [   1      1.5          3      4   ]*(0.324/2) --> Estequiometria
% Reação global: (a+2*k*a)*C3H8 + b'*O2 ----> c*CO2 + (2*c*k)CO + (d+k*d)*H2O + (k*d)*H2
% ========================================================== %

% ======================= Hipóteses ======================== %
% 1) Reação global: (a+2*k*a)*C3H8 + b'*O2 ----> c*CO2 + (2*c*k)CO + (d+k*d)*H2O + (k*d)*H2;
% 2) O N2 não reage no processo;
% 3) Como base de cálculo, temos que a entrada de propano é de 1 mol (a=1);
% ========================================================== %

% ===== Definindo o sistema em forma de matriz: Ax = b ===== %
% Definindo o valor do fator de multiplicação das reações:
k = 0.324/2;

% Definindo a matriz A 4x4 de coeficientes a partir das Equações definidas por balanço:
A = [  1+2*k  ,  0,     0    ,      0;
         0    ,  2,  -(2*k+2),   -(1+k);
     3*(1+2*k),  0,  -(2*k+1),      0;
     8*(1+2*k),  0,     0    ,  -(2+4*k) ];
printf("-> Matriz A de coeficientes:\n"); disp(A);

% Definindo o vetor b 4x1 de solução/termos independentes:
b = [1;
     0;
     0;
     0];
printf("-> Vetor solução b:\n"); disp(b);
% ========================================================== %

% ==================== Método de Cramer ==================== %
function x = metodo_cramer(A, b)                                % Criando a função para cálculo a partir do método de Cramer
  n = length(b);                                                % Variável n com o número de linhas da matriz A
  detA = det(A);                                                % Calculando o determinante da matriz A
  if abs(detA) < 1e-10                                          % Verificando se a matriz apresenta solução única
    disp("A matriz não apresenta solução única");
  end
  x = zeros(n,1);                                               % Cria um vetor x de variáveis a serem substituídas
  for i = 1:n                                                   % Realiza o loop para o cálculo seguindo o método de Cramer
    Ai = A;
    Ai(:,i) = b;
    x(i) = det(Ai)/detA;
  end
end
% ========================================================== %

% ============= Método de Eliminação de Gauss ============== %
function x = metodo_gauss(A, b)
  Ab = [A, b];                                                  % Cria a matriz aumentada
  n = length(b);

  for k = 1:n-1                                                 % Realiza o pivotamento caso necessário
    [~, p] = max(abs(Ab(k:n,k)));                               % Usa a sub coluna do pivô entre as linhas k e n tomando o maior módulo e retornando seu valor máximo
    p = p + k -1;
    if p ~= k                                                   % Realiza a troca das linhas pelo pivô
      temp = Ab(k,:); Ab(k,:) = Ab(p,:); Ab(p,:) = temp;
    end
    for i = k+1:n                                               % Inicia o loop para a eliminação progressiva
      fator = Ab(i,k)/Ab(k,k);
      Ab(i,k:n+1) = Ab(i,k:n+1) - fator*Ab(k,k:n+1);
    end
  end
  x = zeros(n,1);
  for i = n:-1:1                                                % Inicia o loop para a substituição regressiva
    x(i) = (Ab(i,end) - Ab(i,i+1:n)*x(i+1:n))/Ab(i,i);
  end
end
% ========================================================== %

% ================ Método de Matriz Inversa ================ %
function x = metodo_inversa(A, b)
  detA = det(A);
  if abs(detA) < 1e-10
    disp("A matriz não apresenta solução única!");
  end
  A_inv = inv(A);                                               % Calcula a matriz inversa de A
  x = A_inv * b;                                                % Realiza o método de matriz inversa: Ax = b --> A_inv * A * x = A_inv * b --> Ix = A_inv * b --> x = A_inv * b
end
% ========================================================== %

% ================= Chamada das resoluções ================= %
% Solução obtida nos métodos:
x_cramer = metodo_cramer(A, b);
x_gauss = metodo_gauss(A, b);
x_inv = metodo_inversa(A, b);
x_ref = A\b;                                                    % Solução referência do Octave para checar os resultados obtidos

% Vetor solução:
  % xf(1) = C3H8; xf(2) = O2; xf(3) = CO2; xf(4) = CO; xf(5) = H2O; xf(6) = H2
xc = x_cramer;
xcf = [xc(1)+2*xc(1)*k xc(2) xc(3) 2*k*xc(3) xc(4)+k*xc(4) k*xc(4)];

xg = x_gauss;
xgf = [xg(1)+2*xg(1)*k xg(2) xg(3) 2*k*xg(3) xg(4)+k*xg(4) k*xg(4)];

xi = x_inv;
xif = [xi(1)+2*xi(1)*k xi(2) xi(3) 2*k*xi(3) xi(4)+k*xi(4) k*xi(4)];

xr = x_ref;
xrf = [xr(1)+2*xr(1)*k xr(2) xr(3) 2*k*xr(3) xr(4)+k*xr(4) k*xr(4)];

% Composição final:
nt_prop = xrf(1);
nitrog_i = 3.7619*xrf(2);
nitrog_f = 3.7619*xrf(2);
nt_ar = xrf(2) + nitrog_i;                                      % Calculando o número total de mols de reagentes
nt_p = xrf(3) + xrf(4) + xrf(5) + xrf(6) + nitrog_f;            % Calculando o número total de mols de produtos

  % Dividindo a estequiometria pelo número total de mols para obter a composição molar dos componentes:
xc_comp = [(xcf(1)/nt_prop) (xcf(2)/nt_ar) (nitrog_i/nt_ar) (xcf(4)/nt_p) (xcf(3)/nt_p) (xcf(5)/nt_p) (xcf(6)/nt_p) (nitrog_f/nt_p)];
xg_comp = [(xgf(1)/nt_prop) (xgf(2)/nt_ar) (nitrog_i/nt_ar) (xgf(4)/nt_p) (xgf(3)/nt_p) (xgf(5)/nt_p) (xgf(6)/nt_p) (nitrog_f/nt_p)];
xi_comp = [(xif(1)/nt_prop) (xif(2)/nt_ar) (nitrog_i/nt_ar) (xif(4)/nt_p) (xif(3)/nt_p) (xif(5)/nt_p) (xif(6)/nt_p) (nitrog_f/nt_p)];
xr_comp = [(xrf(1)/nt_prop) (xrf(2)/nt_ar) (nitrog_i/nt_ar) (xrf(4)/nt_p) (xrf(3)/nt_p) (xrf(5)/nt_p) (xrf(6)/nt_p) (nitrog_f/nt_p)];

%Print das resoluções:
printf("\nComposição a partir de resolução por Método de Cramer:\n"); disp(xc_comp);
printf("\nComposição a partir de resolução por Método de Eliminação de Gauss:\n"); disp(xg_comp);
printf("\nComposição a partir de resolução por Método de Matriz Inversa:\n"); disp(xi_comp);
printf("\nComposição a partir de resolução por Solução de referência do Octave:\n"); disp(xr_comp);
% ========================================================== %

% =============== Gráficos de Entrada/Saída ================ %
componentes_entrada = {'C3H8', 'O2', 'N2 Entrada'};
componentes_saida = {'CO', 'CO2', 'H2O', 'H2' , 'N2 Saída'};
xr_entrada = xr_comp(1:3);                                     % Valores de entrada - usando a composição de referência
xr_saida = xr_comp(4:8);                                       % Valores de saída - usando a composição de referência
figure;

% Subplot 1 - Entrada de Propano (C3H8):
subplot(2,2,1);
bar(xr_entrada(1), 'FaceColor', [0.2 0.6 0.8]);
title('Composição da Corrente de Entrada');
ylabel('Fração Molar');
set(gca, 'XTickLabel', componentes_entrada(1));
xtickangle(45);
grid on;
ylim([0 1.2]);
text(1, xr_entrada(1)+0.06, sprintf('%.1f%%', xr_entrada(1)*100), 'HorizontalAlignment','center', 'FontSize', 11);

% Subplot 2 - Entrada de Ar seco (O2 e N2):
subplot(2,2,2);
bar(xr_entrada(2:3), 'FaceColor', [0.2 0.6 0.8]);
title('Composição da Corrente de Entrada');
ylabel('Fração Molar');
set(gca, 'XTickLabel', componentes_entrada(2:3));
xtickangle(45);
grid on;
ylim([0 1]);
for i = 1:length(xr_entrada(2:3))
    text(i, xr_entrada(i+1)+0.06, sprintf('%.1f%%', xr_entrada(i+1)*100), 'HorizontalAlignment','center', 'FontSize', 11);
end

% Subplot 3 - Saída:
subplot(2,2,3:4);
bar(xr_saida, 'FaceColor', [0.8 0.4 0.2]);
title('Composição da Corrente de Saída');
ylabel('Fração Molar');
set(gca, 'XTickLabel', componentes_saida);
xtickangle(45);
grid on;
ylim([0 0.85]);
for i = 1:length(xr_saida)
    text(i, xr_saida(i)+0.055, sprintf('%.1f%%', xr_saida(i)*100),'HorizontalAlignment','center', 'FontSize', 11);
end
% ========================================================== %

