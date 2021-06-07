clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simulação
tempo_total = 10*60;   % min
Ts = 1;              % min
nsim = ceil(tempo_total/Ts);
t = Ts : Ts*nsim;

% Estado Estacionário Exato
uss = [156.080602324805,574.692626129846];
y = [0.0965829791616299,3.29484213002587,323.560569143955,303.319300022361,0.000398359512973287,19.2031941791062];

% Obtenção do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','final'));

%% Malha SISO
% Parâmetros de sintonia
PID_Qi = [-83.148611713950790, -12.512685662730720, -0.269884930099017];
PID_Qc = [-2.128358892805603e+02, -16.501882846447092, -4.281661430481993e-04];


%% Resultados
Erro_de_medicao_em_T = 1;
[Y,U,e,setpoint] = malhas_paralelas(reator, yss, uss, uss, nsim, Ts, PID_Qc, PID_Qi, Erro_de_medicao_em_T);

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

kp = reator.Ap*exp(-reator.Ep./T);


Qi = U(:,1);
Qc = U(:,2);

%% Gráficos
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

%% Gráficos Apresentados
figure
subplot(211)
plot(t/60,T,'LineWidth',2)
hold on
stairs(t/60,setpoint(:,2),'r--','LineWidth',2)
if Erro_de_medicao_em_T > 0
    stairs(t/60,setpoint(:,2)+1,'k--','LineWidth',1.5)
    stairs(t/60,setpoint(:,2)-1,'k--','LineWidth',1.5)
end
title('Temperatura do Reator')
legend('Temperatura','Set-point')
ylabel('T/(K)')
xlabel('t/(h)')
xticks([0:1:max(t/60)])
grid on

subplot(212)
plot(t/60,visc,'LineWidth',2)
title('Viscosidade')
ylabel('\eta /(u.v.)')
hold on
stairs(t/60,setpoint(:,1),'r--','LineWidth',2)
legend('Viscosidade','Set-point')
xlabel('t/(h)')
xticks([0:1:max(t/60)])
grid on

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

function [Y,U,e,setpoint] = malhas_paralelas(reator, y0, u0, uss, nsim, Ts, PID_Qc, PID_Qi,Erro_de_medicao_em_T)
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
        if i < 1*60
            % Mantem estacionário
        elseif i < 2*60
            Tsp = 323.56;
            %Tsp = 320;
        elseif i < 5*60
            Viscosp = 2.964;
            Tsp = 323.56;
        elseif i < 7*60 
            Viscosp = 3.17;
            Tsp = 320.8;
        else
           reator.Qm =  302.4;  % 80% do valor inicial
        end
        setpoint = [setpoint; [Viscosp, Tsp]];
        
        
        % Coleta de CV
        e(i+1,1) = Viscosp - reator.viscosidade;    % Viscosidade
        e(i+1,2) = Tsp - y0(3)+ + Erro_de_medicao_em_T*2*(rand()-.5);    % Temperatura
        
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