clc; clear all; close all;

% ========================================================== %
% Objetivo: Resolver o sistema de equações não-lineares pelo método de Newton-Raphson e determinar as concentrações dos reagentes e produtos:
% ========================================================== %

% ======================= Hipóteses ======================== %
% 1) Sistema em estado estacionário:
%   1.1) Sem variação com o tempo (dt = 0);
%   1.2) A vazão de saída e entrada serão iguais (F0 = Ff = F).
% 2) Não há presença do produto C na entrada do reator (Cc0 = 0);
% 3) Reação: A + B --> C com uma constante cinética k.
% ========================================================== %

% ======================== Início ========================== %
% Dados iniciais:
Ca0 = 200;            % mol/m3
Cb0 = 200;            % mol/m3
Cc0 = 0;              % mol/m3
V = 40;               % m3
k = 0.0045863;        % m3/mol.min
F = 3;                % m3/min

% Condições iniciais:
x0 = [Ca0; Cb0; Cc0];                                           % Chute inicial
% ========================================================== %

% =============== Método de Newton-Raphson ================= %
function [x, n_iter, historico] = met_newton_raphson (x, F, V, k, Ca0, Cb0, Cc0)
  tol = 1e-6;                                                   % Definindo a tolerância
  erro = 100;
  max_iter = 150;                                               % Definindo o número máximo de iterações
  n_iter = 0;                                                   % Definindo a partida da iteração

  while n_iter < max_iter && erro > tol                         % Laço para resolução
    % Chute inicial:
    Ca = x(1); Cb = x(2); Cc = x(3);

    % Leis de velocidade:
    ra = -k*Ca*Cb;
    rb = -k*Ca*Cb;
    rc = k*Ca*Cb;

    % Funções:
    fa = F*Ca0 - F*Ca + ra*V;
    fb = F*Cb0 - F*Cb + rb*V;
    fc = F*Cc0 - F*Cc + rc*V;
    f = -[fa; fb; fc];

    % Derivadas das funções:
    dfadCa = -F-k*Cb*V;
    dfadCb = -k*Ca*V;
    dfadCc = 0;
    dfbdCa = -k*Cb*V;
    dfbdCb = -F-k*Ca*V;
    dfbdCc = 0;
    dfcdCa = k*Cb*V;
    dfcdCb = k*Ca*V;
    dfcdCc = -F;

    % Matriz Jacobiana:
    J = [ dfadCa , dfadCb , dfadCc;
          dfbdCa , dfbdCb , dfbdCc;
          dfcdCa , dfcdCb , dfcdCc ];

    % Resolvendo o sistema:
    dx = inv(J)*f;                                             % Resolvendo o sistema linear
    x_new = x + dx;                                            % Atualiza a solução
    erro = max(abs(dx./x));                                    % Calcula o erro para comparar com a tolerância
    x = x_new;                                                 % Atualiza os valores da próxima iteração
    n_iter = n_iter + 1;                                       % Atualiza qual iteração está acontecendo

    % Verificando a convergência:
    historico(n_iter+1,:) = [n_iter, x(1), x(2), x(3), erro];  % Salva os valores de cada iteração
    if erro < tol                                              % Compara o erro calculado com a tolerância definida
      historico = historico(1:n_iter+1,:);                     % Salva qual iteração está acontecendo
      break;
    endif
  endwhile
  if n_iter == max_iter && erro > tol                          % Caso o número máximo de iterações seja atingido
    printf('Atingiu o número máximo de iterações sem atingir a convergência\n')
  endif
endfunction
% ========================================================== %

% ================== Respota do sistema ==================== %
% Definindo as variáveis resposta calculadas pelo método:
[x_final, n_iter_final, historico] = met_newton_raphson(x0, F, V, k, Ca0, Cb0, Cc0);

% Print no command dos resultados finais de iteração e concentrações:
fprintf('- Número de iterações: %d\n', n_iter_final);
fprintf('\n- Concentração inicial de A (Ca0): %.2f mol/m³', x0(1));
fprintf('\n- Concentração final de A (Caf): %.2f mol/m³\n', x_final(1));
fprintf('\n- Concentração inicial de B (Cb0): %.2f mol/m³', x0(2));
fprintf('\n- Concentração final de B (Cbf): %.2f mol/m³\n', x_final(2));
fprintf('\n- Concentração inicial de C (Cc0): %.2f mol/m³', x0(3));
fprintf('\n- Concentração final de C (Ccf): %.2f mol/m³\n', x_final(3));
% ========================================================== %

% ==== Análise de sensibilidade: Variação de Ca0 e Cb0 ===== %
Ca0_sens = linspace(50, 350, 100);
Cb0_sens = linspace(50, 350, 100);
Cc0_sens = zeros(length(Ca0_sens), length(Cb0_sens));
[Ca0_mesh,Cb0_mesh] = meshgrid(Ca0_sens, Cb0_sens);

for i = 1:length(Ca0_sens)
    for j = 1:length(Cb0_sens)
      x0_sens = [Ca0_sens(i), Cb0_sens(j), Cc0];
      [x_sens, n, h] = met_newton_raphson(x0_sens, F, V, k, Ca0_sens(i), Cb0_sens(j), Cc0);
      Cc0_sens(i,j) = x_sens(3);
    endfor
endfor
% ========================================================== %

% ======================== Gráficos ======================== %
componentes = {'A', 'B', 'C'};                                 % Define quais os componentes trabalhados

% ======= Gráfico: Concentrações de Entrada e Saída ======== %
% Agrupando os valores iniciais e finais de concentração por componente:
inicial = [x0(1), x0(2), x0(3)];
final = [x_final(1), x_final(2), x_final(3)];

figure(1);

% Estruturando o gráfico e as legendas:
bar([inicial; final]');
set(gca, 'XTickLabel', componentes);
ylabel('Concentração (mol/m³)');
title('Concentrações de entrada e saída');
legend('[ ] Entrada', '[ ] Saída', 'Location', 'northeast');
ylim([0 210]);

% Apresentando os valores de cada coluna:
  % Coluna de Ca:
text(1-0.2,inicial(1)-25,sprintf('%.1f %s',inicial(1),'mol/m³'),'HorizontalAlignment','center',...
     'VerticalAlignment','middle','FontSize',11,'Rotation',90);
text(1+0.17,final(1)-25,sprintf('%.1f %s',final(1),'mol/m³'),'HorizontalAlignment','center',...
     'VerticalAlignment','middle','FontSize',11,'Rotation',90);

  % Coluna de Cb:
text(2-0.2,inicial(2)-25,sprintf('%.1f %s',inicial(2),'mol/m³'),'HorizontalAlignment','center',...
     'VerticalAlignment','middle','FontSize',11,'Rotation',90);
text(2+0.17,final(2)-25,sprintf('%.1f %s',final(2),'mol/m³'),'HorizontalAlignment','center',...
     'VerticalAlignment','middle','FontSize',11,'Rotation',90);

  % Coluna de Cc:
text(3+0.17,final(3)-25,sprintf('%.1f %s',final(3),'mol/m³'),'HorizontalAlignment','center',...
     'VerticalAlignment','middle','FontSize',11,'Rotation',90);
% ========================================================== %

% ============= Gráfico: Conversão de A e B ================ %
% Calculando a conversão dos reagentes A e B:
conversao_A = (x0(1) - x_final(1)) / x0(1) * 100;
conversao_B = (x0(2) - x_final(2)) / x0(2) * 100;

figure(2);

% Definindo a cor de cada coluna:
bar(1, conversao_A, 'facecolor', [0.2, 0.6, 0.8]);             % Define a cor da coluna A
hold on;
bar(3, conversao_B, 'facecolor', [0.8, 0.2, 0.2]);             % Define a cor da coluna B
hold off;

% Estruturando o gráfico e as legendas:
set(gca, 'XTick', [1, 3]);
set(gca, 'XTickLabel', {'A', 'B'});
ylabel('Conversão (%)');
title('Conversão de A e B');
ylim([0 100]);

% Apresentando os valores de cada coluna:
text(1,conversao_A+3,sprintf('%.1f%%',conversao_A),'HorizontalAlignment','center','FontSize', 11);
text(3,conversao_B+3,sprintf('%.1f%%',conversao_B),'HorizontalAlignment','center','FontSize', 11);
% ========================================================== %

% ========= Gráfico: Concentração por iteração ============= %
figure(3);

% Atualizando os valores iniciais de concentração na iteração inicial (zero):
historico_corrigido = historico;
historico_corrigido(1,2) = Ca0;
historico_corrigido(1,3) = Cb0;
historico_corrigido(1,4) = Cc0;

% Puxando os históricos de concentração a cada iteração e plotando em gráfico:
plot(historico_corrigido(:,1),historico_corrigido(:,2),'ro-','LineWidth', 2,'MarkerSize', 6);
hold on;
plot(historico_corrigido(:,1),historico_corrigido(:,3),'bs-','LineWidth', 2,'MarkerSize', 10);
plot(historico_corrigido(:,1),historico_corrigido(:,4),'g^-','LineWidth', 2,'MarkerSize', 6);
hold off;

% Estruturando o gráfico e as legendas:
xlabel('Número de Iteração');
ylabel('Concentração (mol/m³)');
title('Concentrações ao Longo das Iterações');
legend('Ca', 'Cb', 'Cc', 'Location', 'northeast');
grid on;
% ========================================================== %

% ========== Gráfico: Análise de sensibilidade ============= %
figure(4);

% Subplot 1 - Gráfico 2D da Análise de sensibilidade:
subplot(2,1,1);

  % Puxando os valores da análise de sensibilidade:
contourf(Ca0_mesh, Cb0_mesh, Cc0_sens);

  % Estruturando o gráfico e as legendas:
c = colorbar;
ylabel(c, 'C_C (mol/m³)', 'FontSize', 12);
xlabel('C_{A0} (mol/m³)', 'FontSize', 12);
ylabel('C_{B0} (mol/m³)', 'FontSize', 12);
title('Análise de sensibilidade: Variação da corrente de entrada (2D)');
colormap(pink);
hold on;

  % Plotando o ponto de controle (usa Ca0, Cb0 e Cc definidos no método original):
plot([Ca0,Ca0],[0,Cb0], '--k;;', 'LineWidth', 2);
plot([0,Ca0],[Cb0,Cb0], '--k;;', 'LineWidth', 2);
plot(Ca0, Cb0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
text(Ca0+5, Cb0+5, ['C_C = ', num2str(x_final(3), 4), ' mol/m³'], 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 12);
hold off;

% Subplot 2 - Gráfico 3D da Análise de sensibilidade
subplot(2,1,2);

  % Puxando os valores da análise de sensibilidade:
surf(Ca0_mesh, Cb0_mesh, Cc0_sens, 'EdgeColor', 'none');
shading interp;

  % Estruturando o gráfico e as legendas:
c = colorbar;
ylabel(c, 'C_C (mol/m³)', 'FontSize', 12);
xlabel('C_{A0} (mol/m³)');
ylabel('C_{B0} (mol/m³)');
zlabel('C_C (mol/m³)');
title('Análise de sensibilidade: Variação da corrente de entrada (3D)', 'FontSize', 11, 'FontWeight', 'bold');
colormap(copper);
view(45, 35);
axis tight;
set(gcf, 'Position', [100, 100, 800, 850]);
hold on;

  % Plotando o ponto de controle (usa Ca0, Cb0 e Cc definidos no método original):
plot3(Ca0, Cb0, x_final(3), 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
altura_maxima = max(Cc0_sens(:));
x_offset = 15;                                                 % Deslocamento no eixo x
y_offset = 15;                                                 % Deslocamento no eixo y
z_offset = altura_maxima * 0.20;                               % Deslocamento no eixo z
text(Ca0+x_offset, Cb0+y_offset, x_final(3)+z_offset, ['C_C = ', num2str(x_final(3), 4), ' mol/m³'],'Color', 'white','FontWeight', 'bold',...
     'FontSize', 12,'BackgroundColor', 'black','HorizontalAlignment', 'center','EdgeColor', 'white','Margin', 2);
hold off;
% ========================================================== %

