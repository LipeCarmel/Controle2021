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

% Chute inicial como um controlador P de ganho unit�rio
PID = [1, 0, 0];

