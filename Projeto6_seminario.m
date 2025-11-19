clc; clear all; close all;
% ======================= Seminário ======================== %
% Alunos:
  % Guilherme Henrique da Silva Baltazar - RA: 298888
  % Leonardo Aparecido Ferreira Souza - RA: 268932
% ========================================================== %

% ======================= Objetivo ========================= %
% Resolver o sistema de EDOs proposto no Exemplo 8.7 do livro de "Elemento de Engenharia das Reações
% Químicas" de H. Scott Fogler (6th edição).
% ========================================================== %

% ======================= Contexto ========================= %
% Deseja-se operar um reator semibatelada com duas reações em fase líquida, descritas a seguir:
% Reação 1: A + 2B --> C
% Reação 2: 2A + 3C --> D
% São fornecidos dados iniciais para o processo.
% ========================================================== %

% ======================== Início ========================== %
% Parâmetros do método RK4:
tol = 1e-6;               % Tolerância
max_iter = 1e6;           % Máximo de iterações para o método
h = 0.005;                % Passo (min)

% Parâmetros do processo:
FA0 = 3;                            % mol/min
vazao_inicial = 10;                 % L/min
Vol_min = 1000; Vol_max = 2000;     % L
CB0 = 0.2;                          % mol/L
k1a = 10;                           % [(L^3/mol)^2]/min (Para o componente A)
k2c = 15;                           % [(L^3/mol)^4]/min (Para o componente C)

% Condições inicias:
NA0 = 0; NB0 = 200; NC0 = 0; ND0 = 0;   % mol
N0 = [NA0; NB0; NC0; ND0];              % Vetor de condição inicial
% ========================================================== %

% ================= Reator Semibatelada ==================== %
function dNdt = reator_semib(t, N)
  % Parâmetros do processo:
  FA0 = 3;                            % mol/min
  vazao_inicial = 10;                 % L/min
  Vol_min = 1000; Vol_max = 2000;     % L
  CB0 = 0.2;                          % mol/L
  k1a = 10;                           % [(L^3/mol)^2]/min (Para o componente A)
  k2c = 15;                           % [(L^3/mol)^4]/min (Para o componente C)

  % Calculando o volume do reator ao longo do tempo:
  V = Vol_min + vazao_inicial*t;

  % Condições iniciais:
  NA = N(1); NB = N(2); NC = N(3); ND = N(4);     % mol

  % Calculando as concentrações ao longo do tempo:
  CA = NA/V; CB = NB/V; CC = NC/V; CD = ND/V;     % mol/L

  % Avaliando as leis de taxa de reação:
    % A + 2B --> C (Reação 1):
  r1 = k1a * CA * CB^2;

    % 2A + 3C --> D (Reação 2):
  r2 = k2c * CA^2 * CC^3;

  % Taxas de consumo e formação:
    % Reação 1:
  rA1 = -r1;
  rB1 = -2*r1;
  rC1 = r1;

    % Reação 2:
  rA2 = (-2/3)*r2;
  rC2 = -r2;
  rD2 = (1/3)*r2;

    % Taxas resultantes:
  rA = rA1 + rA2;
  rB = rB1;
  rC = rC1 + rC2;
  rD = rD2;

  % EDOs dos balanços molares:
  dNAdt = rA*V + FA0;
  dNBdt = rB*V;
  dNCdt = rC*V;
  dNDdt = rD*V;

  % Retornando a solução do modelo:
  dNdt = [dNAdt; dNBdt; dNCdt; dNDdt];
endfunction
% ========================================================== %

% ==================== Runge-Kutta 4 ======================= %
function [t, N, n_iter] = mRK4(dNdt, N0, h, tol, max_iter)
  % Condições iniciais:
  N(:,1) = N0; t(1) = 0; n_iter = 0;

  % Restrição de volume máximo:
  Vol_min = 1000; Vol_max = 2000;
  vazao_inicial = 10;

  % Método RK4:
  for i = 1:max_iter
    % Calculando as inclinações:
    k1 = dNdt(t(i), N(:,i));
    k2 = dNdt(t(i)+(h/2), N(:,i)+(h/2)*k1);
    k3 = dNdt(t(i)+(h/2), N(:,i)+(h/2)*k2);
    k4 = dNdt(t(i)+h, N(:,i)+h*k3);

    % Calculando a previsão final no fim do intervalo:
    N_rk = N(:,i) + (h/6) * (k1+(2*k2)+(2*k3)+k4);

    % Aplica a restrição de número de mols não serem negativos:
    N_rk = max(N_rk, 0);

    % Verificação se o volume máximo foi atingido:
    Vol_atual = Vol_min + vazao_inicial*(t(i) + h);
    if Vol_atual >= Vol_max
      N(:,i+1) = N(:,i);
      t(i+1) = t(i);
      n_iter = i;
      break;
    endif

    % Calculando os próximos valores de iteração:
    N(:,i+1) = N_rk; t(i+1) = t(i) + h;

    % Avalia a convergência:
    if i > 1
      erro = norm(N_rk - N(:,i))/norm(N(:,i));
      if erro < tol
        N = N(:,1:i+1); t = t(1:i+1); n_iter = i;
        return;
      endif
    endif
  endfor

  % Aviso caso o número máximo de iterações tenha sido atinjido:
  if n_iter == max_iter
    fprintf('Atingiu o máximo de iterações no método de RK4\n');
  endif
endfunction
% ========================================================== %

% =================== Respota do modelo ==================== %
% Respostas finais geradas no modelo:
[t_RK4, N_RK4, ni_RK4] = mRK4(@reator_semib, N0, h, tol, max_iter);

% Número de mols finais:
NA = N_RK4(1,:); NB = N_RK4(2,:); NC = N_RK4(3,:); ND = N_RK4(4,:);
% ========================================================== %

% =================== Parâmetros finais ==================== %
% Calculando o volume ao longo do tempo:
V = Vol_min + vazao_inicial * t_RK4;

% Calculando a seletividade entre os produtos C e D:
  % Criando a matriz de seletividade:
SeletCD = zeros(size(t_RK4));

  % Calculando os valores de seletividade:
for i = 1:length(t_RK4)
  if t_RK4(i) > 0.0001 && ND(i) > 0
    SeletCD(i) = NC(i)/ND(i);
  else
    SeletCD(i) = 0;
  endif
endfor

% Calculando as concentrações finais de componentes:
CA_final = NA(end)/V(end);
CB_final = NB(end)/V(end);
CC_final = NC(end)/V(end);
CD_final = ND(end)/V(end);

% Calculando as conversões de A e B:
  % Conversão de A:
NA_alimentado = NA0 + FA0 * t_RK4(end);               % Calculando o número de mols de A (NA) alimentados no processo
XA = (NA_alimentado - NA(end))/NA_alimentado * 100;

  % Conversão de B:
XB = (NB0 - NB(end))/NB0 * 100;
% ========================================================== %

% ==================== Respostas finais ==================== %
% Parâmetros ao final do processo:
fprintf('--> Parâmetros finais do processo:\n');
fprintf('- Tempo final: %.2f min\n', t_RK4(end));
fprintf('- Volume final atingido: %.2f L\n', V(end));
fprintf('- Número de iterações: %d\n', ni_RK4);

% Número de mols de cada componente ao final do processo:
fprintf('\n--> Número de mols de cada componente ao final do processo:\n');
fprintf('- NA final: %.5f mols\n', NA(end));
fprintf('- NB final: %.5f mols\n', NB(end));
fprintf('- NC final: %.5f mols\n', NC(end));
fprintf('- ND final: %.5f mols\n', ND(end));

% Concentração de cada componente ao final do processo:
fprintf('\n--> Concentração de cada componente ao final do processo:\n');
fprintf('- CA final: %.5f mol/L\n', CA_final);
fprintf('- CB final: %.5f mol/L\n', CB_final);
fprintf('- CC final: %.5f mol/L\n', CC_final);
fprintf('- CD final: %.5f mol/L\n', CD_final);

% Conversão de A e B ao final do processo:
fprintf('\n--> Conversão final dos reagentes A e B:\n');
fprintf('Conversão de A (XA): %.2f%%\n', XA);
fprintf('Conversão de B (XB): %.2f%%\n', XB);

% Seletividade de C e D ao final do processo:
fprintf('\n--> Seletividade dos produtos C e D ao final do processo: %.3f\n\n', SeletCD(end));
% ========================================================== %

% ================ Análise de sensibilidade ================ %
% Avaliando o sistema em caso de diferentes alimentações para NA0 e NB0:
NA0_sens = linspace(0, 200, 15);
NB0_sens = linspace(50, 350, 15);

% Formando as matrizes de resultado:
NC_sens = zeros(length(NA0_sens), length(NB0_sens));
ND_sens = zeros(length(NA0_sens), length(NB0_sens));
[NA0_mesh, NB0_mesh] = meshgrid(NA0_sens, NB0_sens);

% Loop para resolução do modelo e RK4 para os diferentes valores de NA0 e NB0:
for i = 1:length(NA0_sens)
    for j = 1:length(NB0_sens)

        % Novas condições iniciais:
        N0_sens = [NA0_sens(i); NB0_sens(j); NC0; ND0];

        % Executando o modelo do reator para solução em RK4:
        [t_sens, N_sens, ~] = mRK4(@reator_semib, N0_sens, h, tol, max_iter);

        % Número de mols de C e D finais:
        NC_sens(i,j) = N_sens(3,end);
        ND_sens(i,j) = N_sens(4,end);
    endfor
endfor
% ========================================================== %

% ==== Gráfico 1: Número de mols (N) ao longo do tempo ===== %
figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Plotando os valores de Ni ao longo do tempo:
plot(t_RK4, NA, 'r-', 'LineWidth', 2, 'DisplayName', 'N_A');
hold on;
plot(t_RK4, NB, 'b-', 'LineWidth', 2, 'DisplayName', 'N_B');
plot(t_RK4, NC, 'g-', 'LineWidth', 2, 'DisplayName', 'N_C');
plot(t_RK4, ND, 'm-', 'LineWidth', 2, 'DisplayName', 'N_D');
hold off;

% Ajustando o gráfico (eixos, legenda e título):
xlabel('Tempo (min)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Número de mols (N)', 'FontSize', 10, 'FontWeight', 'bold');
title('Número de mols (N) ao longo do tempo', 'FontSize', 12, 'FontWeight', 'bold');
legend('N_A', 'N_B', 'N_C', 'N_D', 'Location', 'northeast');
ylim([0, 210]);
grid on;
% ========================================================== %

% ====== Gráfico 2: Seletividade ao longo do tempo ========= %
figure(2);
set(gcf, 'Position', [100, 100, 1200, 600]);

% Subplot 1 - Seletividade C/D ao longo do tempo:
subplot(1,2,1);

  % Plotando a seletividade ao longo do tempo:
plot(t_RK4, SeletCD, 'k-', 'LineWidth', 3);

  % Ajustando o gráfico (eixos, legenda e título):
xlabel('Tempo (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Seletividade (C/D)', 'FontSize', 12, 'FontWeight', 'bold');
title('Seletividade S_{CD} ao longo do tempo', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

  % Apresentando o valor final no gráfico:
hold on;

    % Marca o ponto final no gráfico:
plot(t_RK4(end), SeletCD(end), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', ...
     sprintf('%.3f', SeletCD(end)));

     % Adiciona o valor numérico:
text(t_RK4(end), SeletCD(end)+1e26, sprintf('S_C_D: %.3f', SeletCD(end)), 'FontSize', 10, ...
     'FontWeight', 'bold', 'Color', 'red', 'BackgroundColor', 'white', 'EdgeColor', 'red', ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Margin', 2);
hold off;

% Subplot 2 - Zoom no pico da Seletividade C/D ao longo do tempo:
subplot(1,2,2);

  % Plotando a seletividade ao longo do tempo:
plot(t_RK4, SeletCD, 'k-', 'LineWidth', 3);

  % Ajustando o gráfico (eixos, legenda e título):
xlabel('Tempo (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Seletividade (C/D)', 'FontSize', 12, 'FontWeight', 'bold');
title('Pico da seletividade S_{CD} (Zoom)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, 0.05]);    % Limitando o eixo X para observar o pico ('zoom')
grid on;
% ========================================================== %

% === Gráfico 3: Conversão e Conversão ao longo do tempo === %
% Calculando as conversões ao longo do tempo:
XA_dt = zeros(size(t_RK4));
XB_dt = zeros(size(t_RK4));

for i = 1:length(t_RK4)

    % Conversão de A
    NA_alimentado = NA0 + FA0 * t_RK4(i);
    if NA_alimentado > 0
        XA_dt(i) = (NA_alimentado - NA(i)) / NA_alimentado * 100;
    else
        XA_dt(i) = 0;
    endif

    % Conversão de B
    if NB0 > 0
        XB_dt(i) = (NB0 - NB(i)) / NB0 * 100;
    else
        XB_dt(i) = 0;
    endif
endfor

% Conversões finais:
XA_final = XA_dt(end);
XB_final = XB_dt(end);

% Gerando o gráfico:
figure(3);
set(gcf, 'Position', [100, 100, 1200, 600]);

  % Subplot 1 - Gráfico de barras com as conversões:
subplot(1,2,1);
conversoes_finais = [XA_final, XB_final];
componentes = {'X_A', 'X_B'};
barras = bar(conversoes_finais, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % Adicionando valores no subplot 1:
for i = 1:length(conversoes_finais)
    text(i, conversoes_finais(i) + 1, sprintf('%.2f%%', conversoes_finais(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, ...
         'FontWeight', 'bold', 'Color', 'blue');
endfor

    % Ajustando o gráfico (eixos, legenda e título):
set(gca, 'XTickLabel', componentes, 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Conversão (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Conversões de A e B', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, 100]);

  % Subplot 2 - Conversões ao longo do tempo:

    % Plotando a conversão ao longo do tempo:
subplot(1,2,2);
plot(t_RK4, XA_dt, 'r-', 'LineWidth', 3, 'DisplayName', 'X_A');
hold on;
plot(t_RK4, XB_dt, 'b-', 'LineWidth', 3, 'DisplayName', 'X_B');

    % Ajustando o gráfico (eixos, legenda e título):
xlabel('Tempo (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Conversão (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Conversões de A e B ao longo do tempo', 'FontSize', 14, 'FontWeight', 'bold');
legend('X_A', 'X_B', 'Location', 'northeast');
grid on;
% ========================================================== %

% ========== Gráfico 4: Análise de sensibilidade =========== %
figure(4);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Definindo pontos para determinar o número de mols de C e D para diferentes NA0 e NB0:
analise = [ 0,   200;
            100, 200;
            200, 200;
            100, 50;
            200, 350 ];

  % Definindo cores para cada ponto analisado:
cores = ['r', 'g', 'b', 'm', 'c'];
marcadores = ['o', 's', 'd', '^', 'v'];

% Subplot 1 - Variação de NC final de acordo com NA0 e NB0:
subplot(2,1,1);

  % Avaliando os valores da análise de sensibilidade:
contourf(NA0_mesh, NB0_mesh, NC_sens);

  % Ajustando o gráfico (eixos, legenda e título):
c = colorbar;
ylabel(c, 'N_C Final (mol)', 'FontSize', 12);
xlabel('N_{A0} (mol)', 'FontSize', 12);
ylabel('N_{B0} (mol)', 'FontSize', 12);
title('Produção do Produto C - Análise de Sensibilidade', 'FontSize', 14, 'FontWeight', 'bold');
colormap(copper);
hold on;

% Identificando os pontos de análise:
for k = 1:size(analise, 1)
    NA0_ponto = analise(k, 1);
    NB0_ponto = analise(k, 2);

    % Calculando os valores para NC nos pontos de análise:
    [~, idx_NA] = min(abs(NA0_sens - NA0_ponto));
    [~, idx_NB] = min(abs(NB0_sens - NB0_ponto));
    NC_valor = NC_sens(idx_NA, idx_NB);

    % Gerando o plot do ponto com seu valor:
    plot(NA0_ponto, NB0_ponto, marcadores(k), 'MarkerSize', 10, 'MarkerFaceColor', cores(k), ...
         'MarkerEdgeColor', 'black', 'LineWidth', 2);
    text(NA0_ponto+3.5, NB0_ponto+10, sprintf('N_C=%.1f', NC_valor), 'Color', 'black', ...
         'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', cores(k));
end

% Subplot 2 - Variação de ND final de acordo com NA0 e NB0:
subplot(2,1,2);

  % Avaliando os valores da análise de sensibilidade:
contourf(NA0_mesh, NB0_mesh, ND_sens);

  % Ajustando o gráfico (eixos, legenda e título):
c = colorbar;
ylabel(c, 'N_D Final (mol)', 'FontSize', 12);
xlabel('N_{A0} (mol)', 'FontSize', 12);
ylabel('N_{B0} (mol)', 'FontSize', 12);
title('Produção do Produto D - Análise de Sensibilidade', 'FontSize', 14, 'FontWeight', 'bold');
colormap(copper);
hold on;

% Identificando os pontos de análise:
for k = 1:size(analise, 1)
    NA0_ponto = analise(k, 1);
    NB0_ponto = analise(k, 2);

    % Calculando os valores para ND nos pontos de análise:
    [~, idx_NA] = min(abs(NA0_sens - NA0_ponto));
    [~, idx_NB] = min(abs(NB0_sens - NB0_ponto));
    ND_valor = ND_sens(idx_NA, idx_NB);

    % Gerando o plot do ponto com seu valor:
    plot(NA0_ponto, NB0_ponto, marcadores(k), 'MarkerSize', 10, 'MarkerFaceColor', cores(k), ...
         'MarkerEdgeColor', 'black', 'LineWidth', 2);
    text(NA0_ponto+3.5, NB0_ponto+10, sprintf('N_D=%.1f', ND_valor), 'Color', 'black', ...
         'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', cores(k));
end
% ========================================================== %

