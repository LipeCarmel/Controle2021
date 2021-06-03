clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simulação
tempo_total = 210;   % min
Ts = 1;              % min
nsim = ceil(tempo_total/Ts);
t = Ts : Ts*nsim;

% Ponto de Equilíbrio de Referência
uss = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obtenção do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','off'));

%% Malha SISO
% Sintonia inicial como um controlador P de ganho unitário
PID_Qi = [-83.148611713950790, -12.512685662730720, -0.269884930099017];
PID_Qc = [-1.040120716248880e+03, -16.069082517216078, -2.100918309282853e+03];


%% Resultado
[Y,U,e] = malhas_paralelas(reator, yss, uss, uss, nsim, Ts, PID_Qc, PID_Qi);

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


Qi = U(:,1);

figure
stairs(t,Qi)
ylabel('Qi')

Qc = U(:,2);

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
function J = custo(indice_desempenho, reator, yss, uss, nsim, Ts, PID, CV, MV, setpoint)
    [~,~,e] = malha_SISO(reator, yss, uss, uss, nsim, Ts, PID, CV, MV, setpoint);
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
    if u < 0
        u = 0;
    end
end

function [Y,U,e] = malhas_paralelas(reator, y0, u0, uss, nsim, Ts, PID_Qc, PID_Qi)
    Viscosp = 3.17;
    Tsp = 320.8;
    e = zeros(nsim + 1, 2);
    e(1,1) = Viscosp - reator.viscosidade;  % Viscosidade
    e(1,2) = Tsp - y0(3);  % Temperatura

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
        
        if i < nsim/3
        elseif i < nsim*3/4
            Viscosp = 2.964;
            Tsp = 323.56;
        else
            Viscosp = 2.5361;
            %Tsp = 323.56;
            Tsp = 320.8;
        end
            
        
        % Coleta de CV
        e(i+1,1) = Viscosp - reator.viscosidade;    % Viscosidade
        e(i+1,2) = Tsp - y0(3)+ + (rand()-.5)/10;   % Temperatura
        
        % Qi
        u0(1) = controlador_PID(e(1:i+1,1), Ts, PID_Qi, u0(1), uss(1));
        % Qc
        u0(2) = controlador_PID(e(1:i+1,2), Ts, PID_Qc, u0(2), uss(2));
        
    end
end