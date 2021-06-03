clear all; close all; clc
% Criando o reator
reator = ReatorPolimer;
% Ponto de Equil�brio de Refer�ncia
u = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];
% Variando a entrada
Qc_test = 100 : 10 : 700;

X = [];    % Matriz de pontos de equil�brio
s = [];    % Nota de estabilidade

wb = waitbar(0, 'Calculando...');   % Barra de carregamento
n_test = length(Qc_test);
for i = 1 : n_test
    waitbar(i/n_test, wb)       % Atualizando a barra
    uss = [u(1), Qc_test(i)];   % Alterando a entrada
    
    % Buscando pontos de equil�brio pr�ximos � entrada
    [y,~,flag] = fsolve(@(x) reator.derivadas(0,x, uss), y,...
                optimset('Display','off')); 
            
    % Se um ponto de equil�brio for encontrado
    if flag > 0
       X = [X; [y, Qc_test(i)]];        % Armazena informa��o
       g = ModeloLinear(y, u, reator);  % Linearizando no ponto
       
       % A estabilidade � avaliada pela presen�a de apenas polos com parte 
       % real negativa
       s = [s; all(real(pole(g))<0)];   % s = 1 se est�vel
    end
    
end
close(wb)   % Fecha barra de carregamento

%% Pontos de Equilibrio
Qc = X(:, end);
I = X(:,1);
M = X(:,2);
T = X(:,3);
Tc = X(:,4);
D0 = X(:,5);
D1 = X(:,6);
visc = 0.0012*(D0./D1).^0.71;

figure
plot(Qc, I, '.')
ylabel('Iniciador')
xlabel('Qc')

figure
plot(Qc, M, '.')
ylabel('Mon�mero')
xlabel('Qc')

figure
plot(Qc, T, '.')
ylabel('Temperatura')
xlabel('Qc')

figure
plot(Qc, visc, '.')
ylabel('Viscosidade')
xlabel('Qc')

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