clc;
clear;

#Dados do problema
T =  [-98 -73 -48 -18 -8 2 12 22 27 32 37 52 102 127 152 177 202];
Cp =  [0.1286 0.1308 0.1348 0.1402 0.1414 0.1426 0.1452 0.1467 0.1470 0.1493 0.1509 0.1517 0.1623 0.1707 0.1742 0.1835 0.1877];

#Interpolando os valores de -80 °C e 180 °C através do método linear
Cp80 = interp1(T, Cp, -80, 'linear');
Cp180 = interp1(T, Cp, 180, 'linear');
printf("Para o novo intervalo:\n");
printf("O valor de Cp para T = -80 ºC é %.4f cal/(g.ºC)\n", Cp80);
printf("O valor de Cp para T = 180 ºC é %.4f cal/(g.°C)\n\n", Cp180);

#Dados considerando o intervalo de -80 a 180 °C
Tr = [-80 -73 -48 -18 -8 2 12 22 27 32 37 52 102 127 152 177 180];
Cpr =  [0.1302 0.1308 0.1348 0.1402 0.1414 0.1426 0.1452 0.1467 0.1470 0.1493 0.1509 0.1517 0.1623 0.1707 0.1742 0.1835 0.1840];

#Segmentando as integrais para os métodos do Trapézio, Simpson 1/3 e 3/8
I_trapezio1 = ((Tr(2)-Tr(1))*(Cpr(1)+Cpr(2)))/2;
I_trapezio2 = ((Tr(3)-Tr(2))*(Cpr(2)+Cpr(3)))/2;
I_trapezio3 = ((Tr(4)-Tr(3))*(Cpr(3)+Cpr(4)))/2;
I_simp1_3 = ((Tr(8)-Tr(4))*(Cpr(4)+4*(Cpr(5))+2*(Cpr(6))+4*(Cpr(7))+Cpr(8)))/12;
I_simp3_8 = ((Tr(11)-Tr(8))*(Cpr(8)+3*(Cpr(9))+3*(Cpr(10))+Cpr(11)))/8;
I_trapezio4 = ((Tr(12)-Tr(11))*(Cpr(11)+Cpr(12)))/2;
I_trapezio5 = ((Tr(13)-Tr(12))*(Cpr(12)+Cpr(13)))/2;
I_trapezio6 = ((Tr(14)-Tr(13))*(Cpr(13)+Cpr(14)))/2;
I_simp1_3_2 = ((Tr(16)-Tr(14))*(Cpr(14)+4*(Cpr(15))+(Cpr(16))))/6;
I_trapezio7 = ((Tr(17)-Tr(16))*(Cpr(16)+Cpr(17)))/2;

#Somando as integrais dos diferentes métodos
printf("Para o método utilizando a regra do Trapézio, Simpson 1/3 e Simpson 3/8:\n");
printf("A quantidade de calor estimada é %.4f kJ\n\n", S1 = (I_trapezio1 + I_trapezio2 + I_simp1_3 + I_simp3_8 + I_trapezio3 + I_trapezio4 + I_trapezio5 + I_trapezio6 + I_simp1_3_2 + I_trapezio7)*1000*0.00419);

#Segmentando apenas para a regra dos Trapézios (Simples e Múltipla)
I_trap1 = ((Tr(2)-Tr(1))*(Cpr(1)+Cpr(2)))/2;
I_trap2 = ((Tr(3)-Tr(2))*(Cpr(2)+Cpr(3)))/2;
I_trap3 = ((Tr(4)-Tr(3))*(Cpr(3)+Cpr(4)))/2;
I_trap4 = ((Tr(8)-Tr(4))*(Cpr(4)+2*(Cpr(5))+2*(Cpr(6))+2*(Cpr(7))+Cpr(8)))/8;
I_trap5 = ((Tr(11)-Tr(8))*(Cpr(8)+2*(Cpr(9))+2*(Cpr(10))+Cpr(11)))/6;
I_trap6 = ((Tr(12)-Tr(11))*(Cpr(11)+Cpr(12)))/2;
I_trap7 = ((Tr(13)-Tr(12))*(Cpr(12)+Cpr(13)))/2;
I_trap8 = ((Tr(16)-Tr(13))*(Cpr(13)+2*(Cpr(14))+2*(Cpr(15))+Cpr(16)))/6;
I_trap9 = ((Tr(17)-Tr(16))*(Cpr(16)+Cpr(17)))/2;

#Soma dos trapezios
printf("Para o método utilizando apenas a regra do Trapézio:\n");
printf("A quantidade de calor estimada é %.4f kJ\n\n", S2 = (I_trap1 + I_trap2 + I_trap3 + I_trap4 + I_trap5 + I_trap6 + I_trap7 + I_trap8 + I_trap9)*1000*0.00419);

#Aplicação incluindo o polinomio interpolador encontrado na parte 1
I_tra1 = ((Tr(2)-Tr(1))*(Cpr(1)+Cpr(2)))/2;
I_tra2 = ((Tr(3)-Tr(2))*(Cpr(2)+Cpr(3)))/2;
I_tra3 = ((Tr(4)-Tr(3))*(Cpr(3)+Cpr(4)))/2;
I_pol4 = ((Tr(8)-Tr(4))*(14*Cpr(4)+64*(Cpr(5))+24*(Cpr(6))+64*(Cpr(7))+14*Cpr(8)))/180;
I_sim1_3 = ((Tr(11)-Tr(8))*(Cpr(8)+3*(Cpr(9))+3*(Cpr(10))+Cpr(11)))/8;
I_tra4 = ((Tr(12)-Tr(11))*(Cpr(11)+Cpr(12)))/2;
I_tra5 = ((Tr(13)-Tr(12))*(Cpr(12)+Cpr(13)))/2;
I_sim3_8 = ((Tr(16)-Tr(13))*(Cpr(13)+3*(Cpr(14))+3*(Cpr(15))+Cpr(16)))/8;
I_tra6 = ((Tr(17)-Tr(16))*(Cpr(16)+Cpr(17)))/2;

#Soma
printf("Para o método utilizando as regras sugeridas e o desenvolvido na primeira parte do exercício:\n");
printf("A quantidade de calor estimada é %.4f kJ\n", S3 = (I_tra1 + I_tra2 + I_tra3 + I_pol4 + I_sim1_3 + I_tra4 + I_tra5 + I_sim3_8 + I_tra6)*1000*0.00419);




