clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simulação
tempo_total = 15*60;   % min
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
%PID_Qc = [-1.040120716248880e+03, -16.069082517216078,-2.100918309282853e+03]; % simples
%PID_Qc = [-9.046433675921414e+02, -16.064519647501786, -1.704890573504619e+03]; % duRdu, R = diag(1e-4)
PID_Qc = [-2.128358892805603e+02, -16.501882846447092, -4.281661430481993e-04]; % du(R)du, R = diag(1e-3)

%% Resultado
[Y,U,e,setpoint] = malhas_paralelas(reator, yss, uss, uss, nsim, Ts, PID_Qc, PID_Qi);

I = Y(:,1);
M = Y(:,2);
T = Y(:,3);
Tc = Y(:,4);
D0 = Y(:,5);
D1 = Y(:,6);
visc = 0.0012*(D1./D0).^0.71;
kd = reator.Ad*exp(-reator.Ed./T);
kt = reator.At*exp(-reator.Et./T);

P = (2*reator.fi*kd.*I./kt).^0.5;

kp = reator.Ap*exp(-reator.Ep./T);


Qi = U(:,1);
Qc = U(:,2);

%%
figure
stairs(t/60, Qi)
title('Alimentação de iniciador')
ylabel('Qi/(L/h)')
xlabel('t/(h)')


figure
stairs(t/60,Qc)
title('Alimentação da camisa')
ylabel('Qc/(L/h)')
xlabel('t/(h)')

figure
plot(t/60,I)
title('Concentração de Iniciador')
ylabel('I/(mol/L)')
xlabel('t/(h)')

figure
plot(t/60,M)
title('Concentração de Monômero')
ylabel('M/(mol/L)')
xlabel('t/(h)')

figure
plot(t/60,T,'LineWidth',1.5)
hold on
stairs(t/60,setpoint(:,2),'r--','LineWidth',1.5)
stairs(t/60,setpoint(:,2)+1,'k--','LineWidth',.5)
stairs(t/60,setpoint(:,2)-1,'k--','LineWidth',.5)

title('Temperatura do Reator')
ylabel('T/(K)')
xlabel('t/(h)')

figure
plot(t/60,Tc)
title('Temperatura da Camisa')
ylabel('Tc/(K)')
xlabel('t/(h)')

figure
plot(t/60,visc,'LineWidth',1.5)
ylabel('Viscosidade')
hold on
stairs(t/60,setpoint(:,1),'r--','LineWidth',1.5)
xlabel('t/(h)')

figure
plot(t/60,P)
title('Concentração de Polímero')
ylabel('P/(M)')
xlabel('t/(h)')

figure
plot(t/60,D0)
title('Momento Estatístico de Ordem Zero')
ylabel('D0')
xlabel('t/(h)')

figure
plot(t/60,D1)
title('Momento Estatístico de Primeira Ordem')
ylabel('D1')
xlabel('t/(h)')

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
end

function [Y,U,e,setpoint] = malhas_paralelas(reator, y0, u0, uss, nsim, Ts, PID_Qc, PID_Qi)
    setpoint = [];
    % Iniciando em estado estacionário
    Viscosp = reator.viscosidade;
    Tsp = y0(3);
            
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
        
        % Escolha de cenário de simulação
        if i < 2.5*60
            % Mantem estacionário
        elseif i < 5*60
            %     Viscosp = 3.17;
            %     Tsp = 320.8;
            Viscosp = 2.5361;
            Tsp = 323.56;
        elseif i < 7.5*60
            Viscosp = 2.964;
            Tsp = 323.56;
        elseif i < 10*60 
%             Viscosp = 2.5361;
%             Tsp = 323.56;     % Exemplo mantendo T
            %Tsp = 320.8;
            
            
            Viscosp = 3.17;
            Tsp = 320.8;
        else
           reator.Qm =  302.4;  % 80% do valor inicial
        end
        setpoint = [setpoint; [Viscosp, Tsp]];
        
        
        % Coleta de CV
        e(i+1,1) = Viscosp - reator.viscosidade;    % Viscosidade
        e(i+1,2) = Tsp - y0(3)+ + 2*(rand()-.5);   % Temperatura
        
        % Qi
        u0(1) = controlador_PID(e(1:i+1,1), Ts, PID_Qi, u0(1), uss(1));
        if u0(1) < 50
            u0(1) = 50;
        end
        % Qc
        u0(2) = controlador_PID(e(1:i+1,2), Ts, PID_Qc, u0(2), uss(2));
        if u0(2) < 200
            u0(2) = 200;
        end
        
    end
end