clc;
clear all;
% ========================================================== %
% Objetivo: Resolver o sistema de equações lineares pelo método XXX e definir a composição da corrente de saída:
% ========================================================== %

% ======================= Hipóteses ======================== %
% 1) Com a presença de O2 em excesso, foi considerada a reação de combustão completa:
%   1.1) Não há formação de CO e H2 na corrente de saída.
% 2) Reação global: a*C3H8 + b*O2 + c*N2 --> d*CO2 + e*H2O + f*N2;
% 3) Temos que a razão molar entre CO e CO2 na saída é dada por: yCO/yCO2 = d/e = 0.324 --> d = 0.324*e;
% 4) Temos que a razão molar entre N2 e O2 na entrada é dada por: yN2/yO2 = c/b = 3,762 --> c = 3.762*b
% 5) O N2 não reage no processo (f=c);
% 6) Como base de cálculo, temos que a entrada de propano é de 1 mol (a=1).
% ========================================================== %

% ===== Definindo o sistema em forma de matriz: Ax = b ===== %
% Definindo a matriz A 6x6 de coeficientes a partir das Equações definidas por balanço:
A = [ 0,  1,     1  ,  0;
      0,  0,     0  ,  2;
     -2,  1,     2  ,  1;
      0,  1,  -0.324,  0 ];
printf("-> Matriz A de coeficientes:\n"); disp(A);

% Definindo o vetor b 6x1 de solução/termos independentes:
b = [ 3;
      4;
      0;
      0 ];
printf("\n-> Vetor solução b:\n"); disp(b);
% ========================================================== %

% =============== Método de Decomposição LU ================ %

% ========================================================== %
