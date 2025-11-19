clc;
clear all;
% ========================================================== %
% Objetivo: Resolver o sistema de equações lineares pelos métodos e definir a composição da corrente de saída:
%           1) Método de Cramer;
%           2) Método de Eliminação Gaussiana;
%           3) Método de Matriz inversa.
% ========================================================== %

% ======================= Hipóteses ======================== %
% 1) Reação global: a*C3H8 + b*O2 + c*N2 --> d*CO + e*CO2 + f*H2O + g*H2 + h*N2;
% 2) Temos que a razão molar entre CO e CO2 na saída é dada por: yCO/yCO2 = d/e = 0.324 --> d = 0.324*e;
% 3) Temos que a razão molar entre N2 e O2 na entrada é dada por: yN2/yO2 = c/b = 3,762 --> c = 3.762*b
% 4) O N2 não reage no processo (h=c);
% 5) Como base de cálculo, temos que a entrada de propano é de 1 mol (a=1);
% 6) A formação de H2 é pouco favorável, logo, a quantidade molar presente é próxima de 0 (g~0).
% ========================================================== %

% ===== Definindo o sistema em forma de matriz: Ax = b ===== %
% Definindo a matriz A 4x4 de coeficientes a partir das Equações definidas por balanço:
A = [ 0,  1,     1  ,  0;
      0,  0,     0  ,  1;
     -2,  1,     2  ,  1;
      0,  1,  -0.324,  0 ];
printf("-> Matriz A de coeficientes:\n"); disp(A);

% Definindo o vetor b 4x1 de solução/termos independentes:
b = [ 3;
      4;
      0;
      0 ];
printf("\n-> Vetor solução b:\n"); disp(b);
% ========================================================== %

% ==================== Método de Cramer ==================== %
function x = metodo_cramer(A, b)                                % Criando a função para cálculo a partir do método de Cramer
  n = length(b);                                                % Variável n com o mesmo tamanho do vetor b
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
  x = A_inv * b;                                                % Realiza o método de matriz inversa:Ax = b --> A_inv * A * x = A_inv * b --> Ix = A_inv * b --> x = A_inv * b
end
% ========================================================== %

% ================= Chamada das resoluções ================= %
% Solução obtida nos métodos:
x_cramer = metodo_cramer(A, b);
x_gauss = metodo_gauss(A, b);
x_inv = metodo_inversa(A, b);
x_ref = A\b;                                                    % Solução referência do Octave para checar os resultados obtidos

% Composição final:
nt_prop = 1;                                                    % Número total de mols da corrente de entrada de Propano (F)
nitrog_i = 3.7619*x_ref(1);                                     % Número de mols de nitrogênio na corrente de entrada de Ar Seco (A)
nitrog_f = 3.7619*x_ref(1);                                     % Número de mols de nitrogênio na corrente de saída de Produtos (P)
nt_ar = x_ref(1) + nitrog_i;                                    % Número total de mols da corrente de entrada de Ar Seco (A)
h2 = 0;                                                         % Número de mols de H2 na corrente de saída de Produtos (P)
nt_p = x_ref(2) + x_ref(3) + x_ref(4) + nitrog_f;               % Número total de mols da corrente de saída de Produtos (P)

xc_comp = [(1/nt_prop) (x_cramer(1)/nt_ar) (nitrog_i/nt_ar) (x_cramer(2)/nt_p) (x_cramer(3)/nt_p) (x_cramer(4)/nt_p) (h2/nt_p) (nitrog_f/nt_p)];
xg_comp = [(1/nt_prop) (x_gauss(1)/nt_ar) (nitrog_i/nt_ar) (x_gauss(2)/nt_p) (x_gauss(3)/nt_p) (x_gauss(4)/nt_p) (h2/nt_p) (nitrog_f/nt_p)];
xi_comp = [(1/nt_prop) (x_inv(1)/nt_ar) (nitrog_i/nt_ar) (x_inv(2)/nt_p) (x_inv(3)/nt_p) (x_inv(4)/nt_p) (h2/nt_p) (nitrog_f/nt_p)];
xr_comp = [(1/nt_prop) (x_ref(1)/nt_ar) (nitrog_i/nt_ar) (x_ref(2)/nt_p) (x_ref(3)/nt_p) (x_ref(4)/nt_p) (h2/nt_p) (nitrog_f/nt_p)];

% Print das resoluções:
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
title('Composição da Corrente de Propano');
ylabel('Fração Molar');
set(gca, 'XTickLabel', componentes_entrada(1));
xtickangle(45);
grid on;
ylim([0 1.2]);
text(1, xr_entrada(1)+0.06, sprintf('%.1f%%', xr_entrada(1)*100), 'HorizontalAlignment','center', 'FontSize', 11);

% Subplot 2 - Entrada de Ar seco (O2 e N2):
subplot(2,2,2);
bar(xr_entrada(2:3), 'FaceColor', [0.2 0.6 0.8]);
title('Composição da Corrente de Ar Seco');
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


