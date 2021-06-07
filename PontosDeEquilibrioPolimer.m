clear all; close all; clc
%% Sistema
% Criando o reator
reator = ReatorPolimer;

% Ponto de Equil�brio de Refer�ncia
u = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];

% Variando a entrada
Q_test = [300 : 5 : 700];

%% Buscando pontos de equil�brio
X = [];    % Matriz de pontos de equil�brio
s = [];    % Nota de estabilidade

wb = waitbar(0, 'Calculando...');   % Barra de carregamento
n_test = length(Q_test);
for i = 1 : n_test
    waitbar(i/n_test, wb)       % Atualizando a barra
    uss = [u(1), Q_test(i)];   % Alterando a entrada Qc
    %uss = [Q_test(i), u(2)];   % Alterando a entrada Qi
    
    % Buscando pontos de equil�brio pr�ximos � entrada
    [y,~,flag] = fsolve(@(x) reator.derivadas(0,x, uss), y,...
                optimset('Display','off')); 
            
    % Se um ponto de equil�brio for encontrado
    if flag > 0
       X = [X; [y, Q_test(i)]];        % Armazena informa��o
       g = ModeloLinear(y, u, reator);  % Linearizando no ponto
       
       % A estabilidade � avaliada pela presen�a de apenas polos com parte 
       % real negativa
       s = [s; all(real(pole(g))<0)];   % s = 1 se est�vel
    end
    
end
close(wb)   % Fecha barra de carregamento
if all(s)
    disp('Sempre est�vel.')
end
%% Pontos de Obtidos
close all
Q = X(:, end);
I = X(:,1);
M = X(:,2);
T = X(:,3);
Tc = X(:,4);
D0 = X(:,5);
D1 = X(:,6);
visc = reator.vetor_viscosidade(D0,D1);
kd = reator.Ad*exp(-reator.Ed./T);
kt = reator.At*exp(-reator.Et./T);
P = (2*reator.fi*kd.*I./kt).^0.5;

%% Gr�ficos
figure
plot(Q, I, 'b-','LineWidth',1.5)
title('Concentra��o de Iniciador')
ylabel('I/(mol/L)')
xlabel('Qc')

figure
plot(Q, M, 'b-','LineWidth',1.5)
title('Concentra��o de Mon�mero')
ylabel('M/(mol/L)')
xlabel('Qc')
grid on

figure
plot(Q, T, 'b-','LineWidth',1.5)
title('Temperatura do Reator')
ylabel('T/(K)')
xlabel('Qc/(L/h)')
grid on

figure
plot(Q, Tc, 'b-','LineWidth',1.5)
title('Temperatura da Camisa')
ylabel('T/(K)')
xlabel('Qc/(L/h)')
grid on


figure
plot(Q, visc, 'b-','LineWidth',1.5)
ylabel('Viscosidade')
xlabel('Qc/(L/h)')
grid on


figure
plot(Q, P, 'b-','LineWidth',1.5)
title('Concentra��o de Pol�mero')
ylabel('P/(mol/L)')
xlabel('Qc/(L/h)')


%% Gr�fico apresentado
figure
subplot(211)
plot(Q, visc, 'b-','LineWidth',1.5)
title('Viscosidade')
ylabel('\eta /(u.v.)')
xlabel('Qc/(L/h)')
yticks([2.4 : .2 : 3.5])
grid on

subplot(212)
plot(Q, T, 'b-','LineWidth',1.5)
title('Temperatura do Reator')
ylabel('T/(K)')
xlabel('Qc/(L/h)')
yticks([320:2:330])
grid on

function g = ModeloLinear(xss, u, reator)
% N�mero de estados e entradas
nx = length(xss);
nu = length(u);

% Estados e entradas simb�licas
X = sym('X',[1 nx]).';
U = sym('U',[1 nu]).';

% Equa��o diferencial simb�lica
dXdt = reator.derivadas(0, X, U);

% Lineariza��o
Apos = jacobian(dXdt,X);
Bpos = jacobian(dXdt,U);

% Espa�o de estados cont�nuo
Apos = double(subs(Apos, [X; U], [xss'; u']));
Bpos = double(subs(Bpos, [X; U], [xss'; u']));

% Transformada de Laplace e obten��o de fun��es de transfer�ncia
s = tf('s');
g = (s*eye(size(Apos,1)) - Apos)\(Bpos);
end
