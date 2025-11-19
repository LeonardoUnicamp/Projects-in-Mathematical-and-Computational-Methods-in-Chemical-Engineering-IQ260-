clc; clear all; close all;
% ======================= Projeto 4 ======================== %
% Alunos:
  % Guilherme Henrique da Silva Baltazar - RA: 298888
  % Leonardo Aparecido Ferreira Souza - RA: 268932
% ========================================================== %

% ======================= Objetivo ========================= %
% Resolver o sistema de EDOs pelos métodos de Euler (E), Preditor-Corretor (PC) e Runge-Kutta de 4ª Ordem (RK4)
% ========================================================== %

% ======================== Início ========================== %
% Parâmetros:
tol = 1e-6;               % Tolerância
max_iter = 10000;         % Máximo de iterações para os métodos
h = 0.2;                  % Passo (h)

% Condições inicias:
X0 = 0.14; S0 = 9.99; P0 = 0.13;    % g/L
v0 = [X0; S0; P0];                  % Vetor de condição inicial
% ========================================================== %

% ======================== Modelo ========================== %
function dvdt = bioprocesso(t, v)
  % Dados iniciais para o modelo:
  um = 0.325;             % h-1
  Ks = 10.53; Ki = 105;   % g/L
  Y = 0.183;              % g/g de glicose
  m = 0.12; k2 = 0.12;    % g/g.h
  k1 = 0.36;              % g/g

  % Condições iniciais:
  X = v(1); S = v(2); P = v(3);     % g/L

  % Restrições:
  S = max(S, 0);          % A quantidade de substrato final não pode ser menor que zero (S_final => 0)
  X = max(X, 0);          % A quantidade de biomassa final não pode ser menor que zero (X_final => 0)

  % EDOs:
    % Avalia o cescimento para valores maiores ou iguais a zero de substrato e biomassa:
  if S > 0 && X > 0
    mi = um * (S/(S+Ks)) * (Ki/(Ki+S));   % Taxa de crescimento específico (mi)
  else
    mi = 0;                               % Crescimento nulo caso não tenha substrato ou biomassa (mi=0 se X<=0 e S<=0)
  endif

    % Avalia o consumo e produção para valores maiores ou iguais a zero de substrato e biomassa:
  if S > 0 && X > 0
    dxdt = X * mi;                        % Crescimento de Biomassa (X)
    dsdt = -X * ((1/Y)*mi + m);           % Consumo de Substrato (S)
    dpdt = X * (k1*mi + k2);              % Síntese de Produto (P)

    % Caso os valores de substrato ou biomassa cheguem a zero (estado estacionário):
  else
    dxdt = 0; dsdt = 0; dpdt = 0;
  endif

  % Retornando a solução do modelo:
  dvdt = [dxdt; dsdt; dpdt];
endfunction
% ========================================================== %

% ========================= Euler ========================== %
function [t, v, n_iter] = mEuler(dvdt, v0, h, tol, max_iter)
  % Condições iniciais:
  v(:,1) = v0; t(1) = 0;

  % Método de Euler:
  for i = 1:max_iter
    v_euler = v(:,i) + h * dvdt(t(i), v(:,i));

    % Calculando os próximos valores:
    v_euler = max(v_euler, 0);                % Aplica as restrições
    v(:,i+1) = v_euler; t(i+1) = t(i) + h;

    % Avalia a convergência:
    if i > 1
      erro = norm(v_euler - v(:,i))/norm(v(:,i));
      if erro < tol
        v = v(:,1:i+1); t = t(1:i+1); n_iter = i;
        return;
      endif
    endif
  endfor

  % Aviso caso o número máximo de iterações tenha sido atinjido:
  if n_iter == max_iter
    fprintf('Atingiu o máximo de iterações no método de Euler');
  endif
endfunction
% ========================================================== %

% =================== Preditor-Corretor ==================== %
function [t, v, n_iter] = mPC(dvdt, v0, h, tol, max_iter)
  % Condições iniciais:
  v(:,1) = v0; t(1) = 0;

  % Método Preditor-Corretor:
  for i = 1:max_iter
    % Cálculo do preditor:
    v_pred = v(:,i)+h * dvdt(t(i),v(:,i));

    % Cálculo do corretor:
    v_corretor = v(:,i)+(h/2) * (dvdt(t(i),v(:,i))+dvdt(t(i)+h,v_pred));

    % Calculando os próximos valores:
    v_corretor = max(v_corretor, 0);          % Aplica as restrições
    v(:,i+1) = v_corretor; t(i+1) = t(i) + h;

    % Avalia a convergência:
    if i > 1
      erro = norm(v_corretor - v(:,i))/norm(v(:,i));
      if erro < tol
        v = v(:,1:i+1); t = t(1:i+1); n_iter = i;
        return;
      endif
    endif
  endfor

  % Aviso caso o número máximo de iterações tenha sido atinjido:
  if n_iter == max_iter
    fprintf('Atingiu o máximo de iterações no método do Preditor-Corretor');
  endif
endfunction
% ========================================================== %

% ================= Runge-Kutta 4ª Ordem =================== %
function [t, v, n_iter] = mRK4(dvdt, v0, h, tol, max_iter)
  % Condições iniciais:
  v(:,1) = v0; t(1) = 0;

  % Método RK4:
  for i = 1:max_iter
    % Calculando as inclinações:
    k1 = dvdt(t(i), v(:,i));
    k2 = dvdt(t(i)+(h/2), v(:,i)+(h/2)*k1);
    k3 = dvdt(t(i)+(h/2), v(:,i)+(h/2)*k2);
    k4 = dvdt(t(i)+h, v(:,i)+h*k3);

    % Calculando a previsão final no fim do intervalo:
    v_rk = v(:,i) + (h/6) * (k1+(2*k2)+(2*k3)+k4);

    % Calculando os próximos valores:
    v_rk = max(v_rk, 0);                      % Aplica as restrições
    v(:,i+1) = v_rk; t(i+1) = t(i) + h;

    % Avalia a convergência:
    if i > 1
      erro = norm(v_rk - v(:,i))/norm(v(:,i));
      if erro < tol
        v = v(:,1:i+1); t = t(1:i+1); n_iter = i;
        return;
      endif
    endif
  endfor

  % Aviso caso o número máximo de iterações tenha sido atinjido:
  if n_iter == max_iter
    fprintf('Atingiu o máximo de iterações no método de RK4');
  endif
endfunction
% ========================================================== %

% ================== Respota do sistema ==================== %
fprintf('--- Análise das respostas do sistema ---\n');

% Definindo as variáveis resposta calculadas pelo método:
  % O uso do '@' faz com que a função 'bioprocesso' seja utilizada como 'âncora' para rodar os parâmetros em cada função de método
[t_e, v_e, ni_e] = mEuler(@bioprocesso, v0, h, tol, max_iter);
[t_pc, v_pc, ni_pc] = mPC(@bioprocesso, v0, h, tol, max_iter);
[t_rk, v_rk, ni_rk] = mRK4(@bioprocesso, v0, h, tol, max_iter);

% Método de Euler:
fprintf('--> Método de Euler:\n');
fprintf('- Tempo final: %.2f h\n', t_e(end));           % Mostra o tempo final (quando atinge o estado estacionário)
fprintf('- Número de iterações: %d\n', ni_e);           % Mostra o número total de iterações feitas pelo método
fprintf('- Concentrações finais:\n');
fprintf('  ¬Biomassa (X): %.6f g/L\n', v_e(1,end));     % Mostra a concentração final em estado estacionário de biomassa (X) no sistema
fprintf('  ¬Substrato (S): %.6f g/L\n', v_e(2,end));    % Mostra a concentração final em estado estacionário de substrato (S) no sistema
fprintf('  ¬Produto (P): %.6f g/L\n\n', v_e(3,end));    % Mostra a concentração final em estado estacionário de produto (P) no sistema

% Método do Preditor-Corretor:
fprintf('--> Método Preditor-Corretor:\n');
fprintf('- Tempo final: %.2f h\n', t_pc(end));
fprintf('- Número de iterações: %d\n', ni_pc);
fprintf('- Concentrações finais:\n');
fprintf('  ¬Biomassa (X): %.6f g/L\n', v_pc(1,end));
fprintf('  ¬Substrato (S): %.6f g/L\n', v_pc(2,end));
fprintf('  ¬Produto (P): %.6f g/L\n\n', v_pc(3,end));

% Método de Runge-Kutta de 4ª Ordem:
fprintf('--> Método Runge-Kutta 4ª Ordem:\n');
fprintf('- Tempo final: %.2f h\n', t_rk(end));
fprintf('- Número de iterações: %d\n', ni_rk);
fprintf('- Concentrações finais:\n');
fprintf('  ¬Biomassa (X): %.6f g/L\n', v_rk(1,end));
fprintf('  ¬Substrato (S): %.6f g/L\n', v_rk(2,end));
fprintf('  ¬Produto (P): %.6f g/L\n', v_rk(3,end));
% ========================================================== %

% =================== Erros relativos ====================== %
fprintf('\n--- Análise de erro relativo entre os métodos ---\n');
fprintf('(OBS: Não foi analisado o erro relativo de substrato (S) pois todos atingem o resultado final de 0 g/L)\n');

% Erro relativo entre o método de Euler e o método do Preditor-Corretor:
erro_relatX_EPC = abs((v_pc(1,end) - v_e(1,end))/v_pc(1,end));
erro_relatP_EPC = abs((v_pc(3,end) - v_e(3,end))/v_pc(3,end));

fprintf('\n--> Erro relativo entre o método de Euler e o método do Preditor-Corretor:\n');
fprintf('- Biomassa (X): %.6f%%\n', erro_relatX_EPC*100);
fprintf('- Produto (P): %.6f%%\n', erro_relatP_EPC*100);

% Erro relativo entre o método de Euler e o método de Runge-Kutta 4ª Ordem:
erro_relatX_ERK = abs((v_rk(1,end) - v_e(1,end))/v_rk(1,end));
erro_relatP_ERK = abs((v_rk(3,end) - v_e(3,end))/v_rk(3,end));

fprintf('\n--> Erro relativo entre o método de Euler e o método de RK4:\n');
fprintf('- Biomassa (X): %.6f%%\n', erro_relatX_ERK*100);
fprintf('- Produto (P): %.6f%%\n', erro_relatP_ERK*100);

% Erro relativo entre o método do Preditor-Corretor e o método de Runge-Kutta 4ª Ordem:
erro_relatX_PCRK = abs((v_rk(1,end) - v_pc(1,end))/v_rk(1,end));
erro_relatP_PCRK = abs((v_rk(3,end) - v_pc(3,end))/v_rk(3,end));

fprintf('\n--> Erro relativo entre o método do Preditor-Corretor e o método de RK4:\n');
fprintf('- Biomassa (X): %.6f%%\n', erro_relatX_PCRK*100);
fprintf('- Produto (P): %.6f%%\n', erro_relatP_PCRK*100);
% ========================================================== %

% ================== Ajuste dos gráficos =================== %
% Para melhor visualização do estado estacionário no gráfico, foi adicionado +5 horas no tempo final:
  % (tempo final = tempo onde se atingiu o estado estacionário)
t_extra = 5;                                    % Estende o tempo em 5 h
t_plot_euler = t_e(1):h:t_e(end)+t_extra;       % Usando o tempo final do método de Euler
t_plot_pc = t_pc(1):h:t_pc(end)+t_extra;        % Usando o tempo final do método do Preditor-Corretor
t_plot_rk = t_rk(1):h:t_rk(end)+t_extra;        % Usando o tempo final do método de RK4

% Repetindo os últimos valores para prolongar o estado estacionário:
v_e_plot = [v_e, repmat(v_e(:,end), 1, length(t_plot_euler)-length(t_e))];
v_pc_plot = [v_pc, repmat(v_pc(:,end), 1, length(t_plot_pc)-length(t_pc))];
v_rk_plot = [v_rk, repmat(v_rk(:,end), 1, length(t_plot_rk)-length(t_rk))];
% ========================================================== %

% ========= Gráfico 1: Método de Euler (Questão 1) ========= %
figure(1);

% Subplot 1 - Biomassa por tempo:
subplot(3,1,1);
plot(t_plot_euler, v_e_plot(1,:), 'b-', 'LineWidth', 2);
ylabel('X (g/L)');
xlabel('Tempo (h)');
title('Biomassa vs Tempo');
xlim([0 30]);                                                 % Limita os valores do eixo X
grid on;

% Subplot 2 - Substrato por tempo:
subplot(3,1,2);
plot(t_plot_euler, v_e_plot(2,:), 'b-', 'LineWidth', 2);
ylabel('S (g/L)');
xlabel('Tempo (h)');
title('Substrato vs Tempo');
xlim([0 30]);
grid on;

% Subplot 3 - Produto por tempo:
subplot(3,1,3);
plot(t_plot_euler, v_e_plot(3,:), 'b-', 'LineWidth', 2);
ylabel('P (g/L)');
xlabel('Tempo (h)');
title('Produto vs Tempo');
xlim([0 30]);
grid on;
% ========================================================== %

% ===== Gráfico 2: Comparação dos métodos (Questão 2) ====== %
% Plot 1 - Biomassa por tempo:
figure(2);

  % Subplot 1 - Completo:
subplot(3,1,1);
plot(t_plot_euler, v_e_plot(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(1,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(1,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('X (g/L)');
xlabel('Tempo (h)');
title('Biomassa vs Tempo');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([0 1.6]); xlim([0 30]);                       % Limita os valores dos eixos
grid on;

  % Subplot 2 - Zoom 1 (meio da função):
subplot(3,1,2);
plot(t_plot_euler, v_e_plot(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(1,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(1,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('X (g/L)');
xlabel('Tempo (h)');
title('Zoom 1: Intervalo intermediário');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([0.65 0.8]); xlim([12 13]);
grid on;

  % Subplot 3 - Zoom 2 (Perto do estado estacionário):
subplot(3,1,3);
plot(t_plot_euler, v_e_plot(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(1,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(1,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('X (g/L)');
xlabel('Tempo (h)');
title('Zoom 2: Intervalo próximo ao estado estacionário');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([1.48 1.55]); xlim([23 25.6]);
grid on;

% Plot 2 - Substrato por tempo:
figure(3);

  % Subplot 1 - Completo:
subplot(3,1,1);
plot(t_plot_euler, v_e_plot(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(2,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(2,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('S (g/L)');
xlabel('Tempo (h)');
title('Substrato vs Tempo');
legend('Euler', 'Preditor-Corretor', 'RK4');
xlim([0 30]);
grid on;

  % Subplot 2 - Zoom 1 (meio da função):
subplot(3,1,2);
plot(t_plot_euler, v_e_plot(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(2,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(2,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('S (g/L)');
xlabel('Tempo (h)');
title('Zoom 1: Intervalo intermediário');
legend('Euler', 'Preditor-Corretor', 'RK4');
ylim([1.5 4.5]); xlim([16 20]);
grid on;

  % Subplot 3 - Zoom 2 (Perto do estado estacionário):
subplot(3,1,3);
plot(t_plot_euler, v_e_plot(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(2,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(2,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('S (g/L)');
xlabel('Tempo (h)');
title('Zoom 2: Intervalo próximo ao estado estacionário');
legend('Euler', 'Preditor-Corretor', 'RK4');
ylim([0 0.2]); xlim([24.5 25.6]);
grid on;

% Plot 3 - Produto por tempo:
figure(4);

  % Subplot 1 - Completo:
subplot(3,1,1);
plot(t_plot_euler, v_e_plot(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(3,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(3,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('P (g/L)');
xlabel('Tempo (h)');
title('Produto vs Tempo');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([0 3.2]); xlim([0 30]);
grid on;

  % Subplot 2 - Zoom 1 (meio da função):
subplot(3,1,2);
plot(t_plot_euler, v_e_plot(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(3,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(3,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('P (g/L)');
xlabel('Tempo (h)');
title('Zoom 1: Intervalo intermediário');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([1.2 2.2]); xlim([15 20]);
grid on;

  % Subplot 3 - Zoom 2 (Perto do estado estacionário):
subplot(3,1,3);
plot(t_plot_euler, v_e_plot(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_plot_pc, v_pc_plot(3,:), 'r--', 'LineWidth', 2);
hold on;
plot(t_plot_rk, v_rk_plot(3,:), 'k:', 'LineWidth', 2);
hold on;
ylabel('P (g/L)');
xlabel('Tempo (h)');
title('Zoom 2: Intervalo próximo ao estado estacionário');
legend('Euler', 'Preditor-Corretor', 'RK4', 'Location', 'southeast');
ylim([2.8 3.15]); xlim([24.5 25.6]);
grid on;
% ========================================================== %

% =========== Análise de sensibilidade: Ajustes ============ %
% Para as análises de sensibilidade, foi utilizado o método de RK4.
% Para melhor visualização do estado estacionário no gráfico, foi adicionado +5 horas no tempo final:
  % (tempo final = tempo onde se atingiu o estado estacionário)
t_extra = 5;              % (h)
% ========================================================== %

% == Sensibilidade: Variando o Substrato (S) (Questão 3a) == %
fprintf('\n--- Análise de Sensibilidade 1 - Variação do Substrato inicial (S0) ---\n')

% Definindo os valores de variação:
S0_var = [10, 20, 30];    % g/L
fprintf('(Intervalos analisados para S0: 10.0, 20.0 e 30.0 g/L)\n');

% Estrutura do gráfico e solução:
figure(5);

for i = 1:length(S0_var)
  % Condição inicial nova:
  v0_sens1 = [X0; S0_var(i); P0];

  % Resolvendo por RK4 com os valores novos:
  [t_sens1, v_sens1, ni_sens1] = mRK4(@bioprocesso, v0_sens1, h, tol, max_iter);

  % Adicionando o tempo extrar para melhor visualizar o estado estacionário:
  t_plot_sens1 = t_sens1(1):h:t_sens1(end)+t_extra;
  v_sens1_plot = [v_sens1, repmat(v_sens1(:,end), 1, length(t_plot_sens1)-length(t_sens1))];

  % Subplot 1 - Biomassa:
  subplot(3,1,1);
  plot(t_plot_sens1, v_sens1_plot(1,:), 'LineWidth', 2, 'DisplayName', sprintf('S0 = %.1f g/L', S0_var(i)));
  hold on;
  ylabel('X (g/L)');
  xlabel('Tempo (h)');
  title('Biomassa vs Tempo');
  legend('Location', 'northwest');
  xlim([0 30.2]);
  grid on;

  % Subplot 2 - Substrato:
  subplot(3,1,2);
  plot(t_plot_sens1, v_sens1_plot(2,:), 'LineWidth', 2, 'DisplayName', sprintf('S0 = %.1f g/L', S0_var(i)));
  hold on;
  ylabel('S (g/L)');
  xlabel('Tempo (h)');
  title('Substrato vs Tempo');
  legend();
  xlim([0 30.2]);
  grid on;

  % Subplot 3 - Produto:
  subplot(3,1,3);
  plot(t_plot_sens1, v_sens1_plot(3,:), 'LineWidth', 2, 'DisplayName', sprintf('S0 = %.1f g/L', S0_var(i)));
  hold on;
  ylabel('P (g/L)');
  xlabel('Tempo (h)');
  title('Produto vs Tempo');
  legend('Location', 'northwest');
  xlim([0 30.2]);
  grid on;

  % Apresentando os resultados finais para cada caso:
  fprintf('\n--> Aplicando S0 = %.1f g/L:\n', S0_var(i));
  fprintf('- Tempo final: %.2f h\n', t_sens1(end));
  fprintf('- Número de iterações: %d\n', ni_sens1);
  fprintf('- Concentrações finais:\n');
  fprintf('  ¬Biomassa (X): %.6f g/L\n', v_sens1(1,end));
  fprintf('  ¬Substrato (S): %.6f g/L\n', v_sens1(2,end));
  fprintf('  ¬Produto (P): %.6f g/L\n', v_sens1(3,end));
endfor
% ========================================================== %

% == Sensibilidade: Variando a Biomassa (X) (Questão 3b) === %
fprintf('\n--- Análise de Sensibilidade 2 - Variação da Biomassa inicial (X0) ---\n')

% Definindo os valores de variação:
X0_var = [0.15, 0.5, 1];   % g/L
fprintf('(Intervalos analisados para X0: 0.15, 0.50 e 1.00 g/L)\n');

% Estrutura do gráfico e solução:
figure(6);

for i = 1:length(X0_var)
  % Condição inicial nova:
  v0_sens2 = [X0_var(i); S0; P0];

  % Resolvendo por RK4 com os valores novos:
  [t_sens2, v_sens2, ni_sens2] = mRK4(@bioprocesso, v0_sens2, h, tol, max_iter);

  % Adicionando o tempo extrar para melhor visualizar o estado estacionário:
  t_plot_sens2 = t_sens2(1):h:t_sens2(end)+t_extra;
  v_sens2_plot = [v_sens2, repmat(v_sens2(:,end), 1, length(t_plot_sens2)-length(t_sens2))];

  % Subplot 1 - Biomassa:
  subplot(3,1,1);
  plot(t_plot_sens2, v_sens2_plot(1,:), 'LineWidth', 2, 'DisplayName', sprintf('X0 = %.2f g/L', X0_var(i)));
  hold on;
  ylabel('X (g/L)');
  xlabel('Tempo (h)');
  title('Biomassa vs Tempo');
  legend('Location', 'southeast');
  xlim([0 29]);
  grid on;

  % Subplot 2 - Substrato:
  subplot(3,1,2);
  plot(t_plot_sens2, v_sens2_plot(2,:), 'LineWidth', 2, 'DisplayName', sprintf('X0 = %.2f g/L', X0_var(i)));
  hold on;
  ylabel('S (g/L)');
  xlabel('Tempo (h)');
  title('Substrato vs Tempo');
  legend();
  xlim([0 29]);
  grid on;

  % Subplot 3 - Produto:
  subplot(3,1,3);
  plot(t_plot_sens2, v_sens2_plot(3,:), 'LineWidth', 2, 'DisplayName', sprintf('X0 = %.2f g/L', X0_var(i)));
  hold on;
  ylabel('P (g/L)');
  xlabel('Tempo (h)');
  title('Produto vs Tempo');
  legend('Location', 'southeast');
  xlim([0 29]);
  grid on;

  % Apresentando os resultados finais para cada caso:
  fprintf('\n--> Aplicando X0 = %.2f g/L:\n', X0_var(i));
  fprintf('- Tempo final: %.2f h\n', t_sens2(end));
  fprintf('- Número de iterações: %d\n', ni_sens2);
  fprintf('- Concentrações finais:\n');
  fprintf('  ¬Biomassa (X): %.6f g/L\n', v_sens2(1,end));
  fprintf('  ¬Substrato (S): %.6f g/L\n', v_sens2(2,end));
  fprintf('  ¬Produto (P): %.6f g/L\n', v_sens2(3,end));
endfor
% ========================================================== %

