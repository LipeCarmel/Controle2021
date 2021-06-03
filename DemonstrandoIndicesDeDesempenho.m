clear all; close all; clc;

%% Valores de exemplo
Ts = .1;                % Tempo de amostragem
t = [0 : Ts : 5];       % Vetor tempo
y = 2*(1 - exp(-2*t));  % Sinal
ysp = 2.2*ones(size(t));% Sinal de referência

%% Gráfico do exemplo
figure
plot(t,y,'k-')
hold on
plot(t,ysp,'r--')

%% Avaliação de Desempenho
e = y - ysp;

ISE = indice.ISE(e, Ts);
IAE = indice.IAE(e, Ts);
ITSE = indice.ITSE(e, Ts);
ITAE = indice.ITAE(e, Ts);
resultado = table(ISE, IAE, ITSE, ITAE);
disp(resultado)