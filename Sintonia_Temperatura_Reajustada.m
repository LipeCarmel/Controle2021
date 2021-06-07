clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simula��o
tempo_total = 120;   % min
Ts = 1;              % min
nsim = ceil(tempo_total/Ts);
t = Ts : Ts*nsim;

% Ponto de Equil�brio de Refer�ncia
uss = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obten��o do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','off'));

%% Malha SISO
% Sintonia obtida (ser� igual a sa�da do fmincon salva em novoPID)
PID = [-2.128358892805603e+02, -16.501882846447092, -4.281661430481993e-04];

CV = 3; % T
MV = 2; % Qc
setpoint = 320.8;

%% Sintonia
novoPID = fmincon(@(par) custo(@(e, Ts) indice.ITSE(e, Ts), reator, yss, uss, nsim, Ts, par, CV, MV, setpoint),...
                    PID,eye(3),zeros(3,1),[],[],[],[],[]);

%% Resultado
[Y,U,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, novoPID, CV, MV, setpoint);

I = Y(:,1);
M = Y(:,2);
T = Y(:,3);
Tc = Y(:,4);
D0 = Y(:,5);
D1 = Y(:,6);
visc = reator.vetor_viscosidade(D0,D1);
kd = reator.Ad*exp(-reator.Ed./T);
kt = reator.At*exp(-reator.Et./T);
P = (2*reator.fi*kd.*M./kt).^0.5;

Qc = U(:,2);

figure
stairs(t, Qc, 'Linewidth', 1.5)
ylabel('Qc/(L/h)')

% figure
% plot(t,I)
% ylabel('Iniciador')
% 
% figure
% plot(t,M)
% ylabel('Mon�mero')
% 
figure
plot(t,T)
ylabel('Temperatura')

% figure
% plot(t,Tc)
% ylabel('Temperatura da Camisa')
% 
% figure
% plot(t,visc)
% ylabel('Viscosidade')
% 
% figure
% plot(t,P)
% ylabel('Concentra��o de Pol�mero')

%% Fun��es
function J = custo(indice_desempenho, reator, yss, uss, nsim, Ts, PID, CV, MV, setpoint)
    [~,U,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, PID, CV, MV, setpoint);
    J = indice_desempenho(e, Ts) + sum(diff(U(:,MV)).^2)/1e3;
end
function u = controlador_PID(e, Ts, PID, uk_1, uss)
    % O vetor erro (e) deve conter ao menos duas posi��es.
    % PID � assumido como vetor linha contendo os ganhos proporcional,
    % integral e derivativo. Ou seja:
    % Ke = PID(1);
    % Ki = PID(2);
    % Kd = PID(3);
    % Ti = Ke/Ki;
    % Td = Kd/Ke;
    erros = [e(end); trapz(e)*Ts; (e(end)-e(end-1))/Ts];
    %erros = [e(end); 0; 0];
    
    u = PID*erros + uss;	% em desvio
    duk = u - uk_1;
    if abs(duk) > abs(uss)*.1
       u =  uk_1 + sign(duk)*uss*.1;
    end
end

function [Y,U,e] = malha_SISO(reator, y0, u0, uss, nsim, Ts, PID, CV, MV, setpoint)
    e = zeros(nsim + 1, 1); e(1) = setpoint - y0(CV);
    Y = zeros(nsim, length(y0));
    U = zeros(nsim, length(u0));
    
    for i = 1 : nsim
        % Armazenda valores
        Y(i, :) = y0;
        U(i, :) = u0;
        
        % Simula��o de um instante de amostragem
        [~,y] = ode15s(@(t,y) reator.derivadas(t,y,u0), [0 Ts], y0);
        
        % Atualiza��o
        y0 = y(end, :);
        
        % Coleta de CV
        e(i+1) = setpoint - y0(CV);
        u0(MV) = controlador_PID(e(1:i+1), Ts, PID, u0(MV), uss(MV));
        if u0(MV) < 0
           u0 = 0; 
        end
    end
end