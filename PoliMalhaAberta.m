clear all; close all; clc
% Criando o reator
reator = ReatorPolimer;

% Ponto de equilíbrio
u = [108, 471.6];
y0 = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];

% Alteração de exemplo: degrau em Qi
u = [208, 471.6];
y0 = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];


[t,y] = ode15s(@(t,y) reator.derivadas(t,y,u), [0: .1 : 200], y0);
%%

            %X = [I M T Tc D0 D1]
figure
plot(t,y(:,1))
ylabel('Iniciador')
figure
plot(t,y(:,2))
ylabel('Monômero')
figure
plot(t,y(:,3))
ylabel('Temperatura')
figure
plot(t,y(:,4))
ylabel('Temperatura da Camisa')
figure
plot(t,(y(:,end-1)./y(:,end)).^0.71)
ylabel('Viscosidade')
