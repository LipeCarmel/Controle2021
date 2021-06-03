clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simulação
tempo_total = 80;   % min
Ts = 1;             % min
nsim = ceil(tempo_total/Ts);
t = Ts : Ts*nsim;

% Ponto de Equilíbrio de Referência
uss = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obtenção do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','off'));

%% Malha SISO
% Sintonia inicial como um controlador P de ganho unitário
PID = [-1, 0, 0];
% O vetor contém ganhos proporcional, integral e derivativo. 

CV = 3; % T
MV = 2; % Qc
setpoint = 320.8;
[Y,U,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, PID, CV, MV, setpoint);

%% Gráficos

I = Y(:,1);
M = Y(:,2);
T = Y(:,3);
Tc = Y(:,4);
D0 = Y(:,5);
D1 = Y(:,6);
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
plot(t,(D0./D1).^0.71)
ylabel('Viscosidade')

figure
plot(t,P)
ylabel('Concentração de Polímero')

%% Funções
function u = controlador_PID(e, Ts, PID)
    % O vetor erro (e) deve conter ao menos duas posições.
    % PID é assumido como vetor linha contendo os ganhos proporcional,
    % integral e derivativo. Ou seja:
    % Ke = PID(1);
    % Ki = PID(2);
    % Kd = PID(3);
    % Ti = Ke/Ki;
    % Td = Kd/Ke;
    erros = [e(end); trapz(e)*Ts; (e(end)-e(end-1))/Ts];
    u = PID*erros;  % em desvio
end

function [Y,U,e] = malha_SISO(reator, y0, u0, uss, nsim, Ts, PID, CV, MV, setpoint)
    e = zeros(nsim + 1, 1);
    Y = zeros(nsim, length(y0));
    U = zeros(nsim, length(u0));
    
    for i = 1 : nsim
        % Armazenda valores
        Y(i, :) = y0;
        U(i, :) = u0;
        
        % Simulação de um instante de amostragem
        [~,y] = ode15s(@(t,y) reator.derivadas(t,y,u0), [0 Ts], y0);
        
        % Atualização
        y0 = y(end, :);
        
        % Coleta de CV
        e(i+1) = setpoint - y0(CV);
        u0(MV) = controlador_PID(e(1:i+1), Ts, PID) + uss(MV);
    end
end