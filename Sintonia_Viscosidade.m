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
viscosiss = reator.viscosidade;

%% Malha SISO
% Sintonia
PID = [-83.148611713950790, -12.512685662730720, -0.269884930099017];

MV = 1; % Qi
setpoint = 3.17;

%% Sintonia
novoPID = fmincon(@(par) custo(@(e, Ts) indice.ITAE(e, Ts), reator, yss, uss, nsim, Ts, par, MV, setpoint),...
                    PID, eye(3),zeros(3,1),[],[],[],[],[], optimoptions('fmincon', 'StepTolerance', 1e-7));


%% Resultado
[Y,U,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, novoPID, MV, setpoint);

I = Y(:,1);
M = Y(:,2);
T = Y(:,3);
Tc = Y(:,4);
D0 = Y(:,5);
D1 = Y(:,6);
visc = reator.vetor_viscosidade(D0,D1);
kd = reator.Ad*exp(-reator.Ed./T);
kt = reator.At*exp(-reator.Et./T);
P = (2*reator.fi*kd.*I./kt).^0.5;
Qi = U(:,1);
Qc = U(:,2);

%% Gráficos
figure
stairs(t,Qi)
ylabel('Qi')

figure
stairs(t,Qc)
ylabel('Qc')

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

%% Funções
function J = custo(indice_desempenho, reator, yss, uss, nsim, Ts, PID, MV, setpoint)
    [~,~,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, PID, MV, setpoint);
    J = indice_desempenho(e, Ts);
end
function u = controlador_PID(e, Ts, PID, uk_1, uss)
    % O vetor erro (e) deve conter ao menos duas posições.
    % PID é assumido como vetor linha contendo os ganhos proporcional,
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

function [Y,U,e] = malha_SISO(reator, y0, u0, uss, nsim, Ts, PID, MV, setpoint)
    e = zeros(nsim + 1, 1); e(1) = setpoint - reator.viscosidade;
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
        e(i+1) = setpoint - reator.viscosidade;
        u0(MV) = controlador_PID(e(1:i+1), Ts, PID, u0(MV), uss(MV));
        if u0(MV) < 0
           u0 = 0; 
        end
    end
end