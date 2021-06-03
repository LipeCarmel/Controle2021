clear all; close all; clc

% Criando o reator
reator = ReatorPolimer;

% Tempo de simulação
time = 0 : 1 : 80;


% Ponto de Equilíbrio de Referência
uss = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Obtenção do ponto exato
yss = fsolve(@(x)reator.derivadas(0,x, uss), y, optimset('Display','off'));

% Chute inicial como um controlador P de ganho unitário
PID = [1, 0, 0];

