clc;
clear;

%Definindo um conjunto de pontos aleatórios
%Cálculo do Trabalho Realizado por uma força variável

D = [0 0.5 1 1.5 2 2.5 3]; %Deslocamento [m]
F = [2 1.683 1.327 1.0334 0.847 0.771 0.779]; %Força [N]

%Definindo o passo h:
%Temos que h = (xn - x0)/n.
h = (D(7)-D(1))/6;

%Aplicação do método de integração numérica a partir do polinômio interpolador de grau 4:
T = (h/45) * (14*F(1)+64*F(2)+24*F(3)+64*F(4)+24*F(5)+64*F(6)+14*F(7));

%Aplicação do método de Simpson 1/3 múltipla:
T2 = (h/3) * (F(1)+4*F(2)+2*F(3)+4*F(4)+2*F(5)+4*F(6)+F(7));

%Resultados:
printf("Resultado obtido pelo método integral do polinômio de grau 4: %.4f Joules", T);

printf("\n\nResultado obtido pelo método de Simpson 1/3: %.4f Joules", T2);

