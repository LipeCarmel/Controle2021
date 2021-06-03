clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simula��o
time = 0 : 1 : 80;

% Ponto de Equil�brio de Refer�ncia
uss = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obten��o do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','off'));

%% Sintonia
% Sintonia inicial como um controlador P de ganho unit�rio
PID = [1, 0, 0];
% O vetor cont�m ganhos proporcional, integral e derivativo. 

function u = controlador_PID(e, Ts, PID)
    % O vetor erro (e) deve conter ao menos duas posi��es.
    % PID � assumido como vetor linha contendo os ganhos proporcional,
    % integral e derivativo. Ou seja:
    % Ke = PID(1);
    % Ki = PID(2);
    % Kd = PID(3);
    % Ti = Ke/Ki;
    % Td = Kd/Ke;
    erros = [e(end); trapz(e)*Ts; (e(end)-e(end-1))/Ts];
    u = PID*erros;
end