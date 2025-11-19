clc; clear all; close all;
% ======================= Projeto 5 ======================== %
% Alunos:
  % Guilherme Henrique da Silva Baltazar - RA: 298888
  % Leonardo Aparecido Ferreira Souza - RA: 268932
% ========================================================== %

% ======================= Objetivo ========================= %
% Resolver a EDP da Equação de Laplace e Equação de Calor de Fourier para calcular a distribuição de
% temperatura e fluxo de calor em uma placa de alumínio
% ========================================================== %

% ======================== Início ========================== %
% Parâmetros:
tol = 1e-6;               % Tolerância
max_iter = 10000;         % Máximo de iterações para o método
Lx = 80; Ly = 50;         % cm
k = 0.49;                 % cal/s.cm.ºC

% Criando um menu para escolher o tamanho do passo de discretização:
  % Menu para discretização em x:
dxmenu = menu('Escolha o tamanho do passo de discretização para dx:', ...
       '1 cm', '2 cm', '4 cm', '5 cm', '8 cm', '10 cm', '16 cm', '20 cm', '40 cm', '80 cm');
Dx = [1, 2, 4, 5, 8, 10, 16, 20, 40, 80];

  % Menu para discretização em y:
dymenu = menu('Escolha o tamanho do passo de discretização para dy:', ...
       '1 cm', '2 cm', '5 cm', '10 cm', '25 cm', '50 cm');
Dy = [1, 2, 5, 10, 25, 50];

% Discretização:
dx = Dx(dxmenu);          % cm
dy = Dy(dymenu);          % cm
nx = Lx/dx + 1;
ny = Ly/dy + 1;
T = zeros(ny, nx);        % ºC

% Condições de contorno:
T(1,:) = 100;             % Parte inferior (ºC) (y = 0)
T(end,:) = 60;            % Parte superior (ºC) (y = Ly)
T(:,1) = 25;              % Parte esquerda (ºC) (x = 0)
T(:,end) = 25;            % Parte direita (ºC)  (x = Lx)
% ========================================================== %

% =============== Equação de Gauss-Seidel ================== %
function [T, n_iter] = gauss_seidel(T, tol, max_iter, dx, dy)
  [ny, nx] = size(T);
  n_iter = 0;

  % Iniciando o método de Gauss-Seidel:
  for n_iter = 1:max_iter
    T_anterior = T;

    % Atualiza o valor dos pontos internos a cada iteração:
    for j = 2:ny-1
      for i = 2:nx-1
        T(j,i) = (dy^2 * (T(j, i-1) + T(j, i+1)) + dx^2 * (T(j-1, i) + T(j+1, i)))/(2*dx^2+2*dy^2);
      endfor
    endfor

    % Avaliando a convergência do sistema:
    if n_iter > 1
      erro = max(max(abs(T - T_anterior)));
      if erro < tol;
        return;
      endif
    endif
  endfor

  % Aviso caso o número máximo de iterações tenha sido atingido:
  if n_iter == max_iter
    fprintf('Número máximo de iterações atingido\n');
  endif
endfunction
% ========================================================== %

% ============== Equação de Fluxo (Fourier) ================ %
function [qx, qy, q_result, angulo] = fluxo(T, k, dx, dy)
  [ny, nx] = size(T);
  qx = zeros(ny, nx);
  qy = zeros(ny, nx);

  % Calculando qx:
  for j = 1:ny
    for i = 1:nx

      % Calculando qx para as bordas:
      if i == 1
        qx(j,i) = -k * (T(j, i+1) - T(j, i))/dx;
      elseif i == nx
        qx(j,i) = -k * (T(j, i) - T(j, i-1))/dx;

      % Calculando qx para os pontos internos:
      else
        qx(j,i) = -k * (T(j, i+1) - T(j, i-1))/(2*dx);
      endif
    endfor
  endfor

  % Calculando qy:
  for j = 1:ny
    for i = 1:nx

      % Calculando qy para as bordas:
      if j == 1
        qy(j,i) = -k * (T(j+1, i) - T(j, i))/dy;
      elseif j == ny
        qy(j,i) = -k * (T(j, i) - T(j-1, i))/dy;

      % Calculando qy para os pontos internos:
      else
        qy(j,i) = -k * (T(j+1, i) - T(j-1, i))/(2*dy);
      endif
    endfor
  endfor

  % Calculando o fluxo resultante:
  q_result = sqrt(qx.^2 + qy.^2);

  % Calculando o ângulo:
  angulo = zeros(ny, nx);
  for j = 1:ny
    for i = 1:nx
      if q_result(j,i) > 0

        % Calculando o ângulo caso qx > 0:
        if qx(j,i) > 0
          angulo(j,i) = atan(qy(j,i) / qx(j,i)) * 180/pi;

        % Calculando o ângulo caso qx < 0:
        elseif qx(j,i) < 0
          angulo(j,i) = (atan(qy(j,i) / qx(j,i)) + pi) * 180/pi;

        % Calculando o ângulo para qx = 0:
        else

          % Calculando o ângulo para qx = 0 e qy > 0:
          if qy(j,i) > 0
            angulo(j,i) = 90;

          % Calculando o ângulo para qx = 0 e qy < 0:
          elseif qy(j,i) < 0
            angulo(j,i) = 270;

          % Calculando o ângulo para qx = 0 e qy = 0:
          else
            angulo(j,i) = 0;
          endif
        endif

        % Garantir ângulo entre 0° e 360°:
        if angulo(j,i) < 0
          angulo(j,i) = angulo(j,i) + 360;
        endif
        if angulo(j,i) >= 360
          angulo(j,i) = angulo(j,i) - 360;
        endif
      endif
    endfor
  endfor
endfunction
% ========================================================== %

% ================== Respota do sistema ==================== %
% Solução do método de Gauss-Seidel/Liebmann:
[T_solucao, n_iter] = gauss_seidel(T, tol, max_iter, dx, dy);
fprintf('--> Método de Gauss-Seidel/Liebmann:\n');
fprintf('- Número de iterações: %d\n', n_iter);

% Solução da equação de fluxo de calor (Fourier):
[qx, qy, q_result, angulo] = fluxo(T_solucao, k, dx, dy);
% ========================================================== %

% ===== Gráfico 1: Distribuição do calor (Questão A) ======= %
figure(1);
set(gcf, 'Position', [100, 100, 1200, 500]);

% Subplot 1 - Distribuição do calor 1 (Mapa de calor):
subplot(5,2,[1,3,5,7]);

  % Definindo as dimensões da figura:
x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);
[X, Y] = meshgrid(x, y);

  % Coloração da figura:
contourf(X, Y, T_solucao, 20, 'LineColor', 'none'); colormap(jet);
hold on;
caxis([25 100]);

  % Enquadramento da figura e nome dos eixos:
rectangle('Position', [0, 0, Lx, Ly], 'EdgeColor', 'k', 'LineWidth', 2);
xlabel('Largura (cm)', 'FontSize', 14); ylabel('Altura (cm)', 'FontSize', 14);
set(gca, 'FontSize', 14);
axis equal;

  % Condições de contorno descritas na figura:
text(-10, Ly/2, '25 °C', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue', 'BackgroundColor', 'white');
text(Lx+5, Ly/2, '25 °C', 'Rotation', 270, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue', 'BackgroundColor', 'white');
text(Lx/2, -10, '100 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'red', 'BackgroundColor', 'white');
text(Lx/2, Ly+3, '60 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.85 0.33 0.1], 'BackgroundColor', 'white');
hold off;

% Subplot 2 - Distribuição do calor 2 (Mapa com valores):
subplot(5,2,[2,4,6,8]);

  % Cria o gráfico onde cada elemento da matriz solução é um quadrado:
imagesc(flipud(T_solucao));

  % Coloração da figura:
colormap(jet);
caxis([25 100]);

  % Configuração dos eixos:
xticks(1:nx); yticks(1:ny);
xticklabels(arrayfun(@(i) sprintf('%.0f', (i-1)*dx), 1:nx, 'UniformOutput', false));
yticklabels(arrayfun(@(j) sprintf('%.0f', (j-1)*dy), 1:ny, 'UniformOutput', false));
hold on;

  % Adiciona os valores numéricos em cada quadrado de acordo com a solução:
for j = 1:ny
    for i = 1:nx

        % Define cores para os valores de temperatura:
        if T_solucao(j,i) > 80
            cor = 'white';
        elseif T_solucao(j,i) < 30
            cor = 'white';
        else
            cor = 'black';
        endif
        text(i, ny-j+1, sprintf('%.0f', T_solucao(j,i)), 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', cor);
    endfor
endfor

  % Enquadramento da figura:
rectangle('Position', [0.5, 0.5, nx, ny], 'EdgeColor', 'black', 'LineWidth', 2, 'FaceColor', 'none');
hold off; axis equal; axis off;

% Subplot 3 - Colorbar (para melhor visualização):
subplot(5,2,9:10);
colormap(jet);
caxis([25 100]);
set(gca, 'Visible', 'off');
c = colorbar('southoutside');
set(c, 'XTick', [25:15:100]);
xlabel(c, 'Temperatura (°C)');
set(c, 'fontsize', 14);
pos = get(c, 'position'); pos = [0.32,0.15,0.4,0.03];
set(c, 'position', pos);
% ========================================================== %

% ======== Gráfico 2: Fluxo de calor (Questão B) =========== %
% Distribuição dos fluxos de calor na placa:
figure(2);
set(gcf, 'Position', [100, 100, 1200, 500]);

  % Definindo as dimensões da figura:
x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);
[X, Y] = meshgrid(x, y);

  % Enquadramento da figura:
mask = true(size(X));
mask(:,1) = false; mask(:,end) = false;
rectangle('Position', [0, 0, Lx, Ly], 'FaceColor', 'white', 'EdgeColor', 'k', 'LineWidth', 2);
hold on;

  % Vetores de fluxo:
quiver(X(mask), Y(mask), qx(mask), qy(mask), 1, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);

  % Configurações dos pontos, eixos e do título:
plot(X, Y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
title('Distribuição de Fluxo de Calor na Placa', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Largura (cm)', 'FontSize', 12); ylabel('Altura (cm)', 'FontSize', 12);
axis equal;
xlim([-15, Lx+15]); ylim([-10, Ly+10]);
grid off;

  % Condições de contorno descritas na figura:
text(-3, Ly/2, '25 °C', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(Lx+3, Ly/2, '25 °C', 'Rotation', 270, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(Lx/2, -3, '100 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'red');
text(Lx/2, Ly+3, '60 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.85 0.33 0.1]);
hold off;

% Matrizes de Fluxo de calor e Ângulo em cada ponto:
figure(3);
set(gcf, 'Position', [100, 100, 1200, 800]);

  % Subplot 1 - Matriz do fluxo de calor (Fourier) resultante:
subplot(1,2,1);

    % Enquadramento da figura:
rectangle('Position', [0.5, 0.5, nx, ny], 'FaceColor', 'white', 'EdgeColor', 'none');
hold on;

    % Formatação dos quadrados na figura:
for j = 1:ny
    for i = 1:nx
        rectangle('Position', [i-0.5, j-0.5, 1, 1], 'FaceColor', 'white', ...
                  'EdgeColor', 'black', 'LineWidth', 1);
        text(i, j, sprintf('%.2f', q_result(j,i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'FontWeight', 'bold', 'Color', 'black');
    endfor
endfor

    % Condições de contorno descritas na figura:
text(0, ny/2, '25 °C', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(nx+1, ny/2, '25 °C', 'Rotation', 270, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(nx/2, 0.2, '100 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'red');
text(nx/2, ny+0.8, '60 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.85 0.33 0.1]);

    % Configurações dos eixos e do título:
title('Fluxos de calor resultantes (cal/s·cm²)', 'FontSize', 12, 'FontWeight', 'bold', ...
      'Position', [nx/2, ny+1]);
xticks(1:nx); yticks(1:ny);
xticklabels(arrayfun(@(i) sprintf('%.0f', (i-1)*dx), 1:nx, 'UniformOutput', false));
yticklabels(arrayfun(@(j) sprintf('%.0f', (j-1)*dy), ny:-1:1, 'UniformOutput', false));
axis equal;
xlim([0.5, nx+0.5]); ylim([0.5, ny+0.5]);
grid off; hold off; axis off;

  % Subplot 2 - Matriz do ângulo resultante:
subplot(1,2,2);

    % Enquadramento da figura:
rectangle('Position', [0.5, 0.5, nx, ny], 'FaceColor', 'white', 'EdgeColor', 'none');
hold on;

    % Formatação dos quadrados na figura:
for j = 1:ny
    for i = 1:nx
        rectangle('Position', [i-0.5, j-0.5, 1, 1], 'FaceColor', 'white', ...
                  'EdgeColor', 'black', 'LineWidth', 1);
        text(i, j, sprintf('%.0f', angulo(j,i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'FontWeight', 'bold', 'Color', 'black');
    endfor
endfor

    % Condições de contorno descritas na figura:
text(0, ny/2, '25 °C', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(nx+1, ny/2, '25 °C', 'Rotation', 270, 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'blue');
text(nx/2, 0.2, '100 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'red');
text(nx/2, ny+0.8, '60 °C', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.85 0.33 0.1]);

    % Configurações dos eixos e do título:
title('Ângulos resultantes (graus)', 'FontSize', 12, 'FontWeight', 'bold', ...
      'Position', [nx/2, ny+1]);
xticks(1:nx); yticks(1:ny);
xticklabels(arrayfun(@(i) sprintf('%.0f', (i-1)*dx), 1:nx, 'UniformOutput', false));
yticklabels(arrayfun(@(j) sprintf('%.0f', (j-1)*dy), ny:-1:1, 'UniformOutput', false));
axis equal;
xlim([0.5, nx+0.5]); ylim([0.5, ny+0.5]);
grid off; hold off; axis off;
% ========================================================== %

