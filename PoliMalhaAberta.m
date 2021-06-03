clear all; close all; clc
% Criando o reator
reator = ReatorPolimer;

% Ponto de equilíbrio
u = [108, 471.6];
y0 = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obtenção do ponto exato
y0 = fsolve(@(x)reator.derivadas(0, x, u), y0, optimset('Display','off'));

u = [108, 471.6];
u(1) = 105;


[t,y] = ode15s(@(t,y) reator.derivadas(t,y,u), [0: .1 : 200], y0);

%% y = [I M T Tc D0 D1]

I = y(:,1);
M = y(:,2);
T = y(:,3);
Tc = y(:,4);
D0 = y(:,5);
D1 = y(:,6);

visc = 0.0012*(D1./D0).^0.71;
kd = reator.Ad*exp(-reator.Ed./T);
kt = reator.At*exp(-reator.Et./T);
P = (2*reator.fi*kd.*M./kt).^0.5;

figure
plot(t,I)
ylabel('Iniciador')

figure
plot(t,M)
ylabel('Monômero')

figure
plot(t,T)
ylabel('Temperatura')

figure
plot(t,Tc)
ylabel('Temperatura da Camisa')

figure
plot(t,visc)
ylabel('Viscosidade')

figure
plot(t,P)
ylabel('Concentração de Polímero')


