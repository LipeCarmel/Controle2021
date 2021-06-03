clear all
close all
clc

u = [108, 471.6];
y = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];

Qc_test = 100 : 10 : 700;

X = [];
s = [];

wb = waitbar(0,'Calculando...');
n_test = length(Qc_test);
for i = 1 : n_test
    waitbar(i/n_test, wb)
    uss = [u(1), Qc_test(i)];
    [y,~,flag] = fsolve(@(x) derivadas(0,x, uss), y, optimset('Display','off')); 
    if flag > 0
       X = [X; [y, Qc_test(i)]]; 
       g = ModeloLinear(y,u);
       s = [s; all(real(pole(g))<0)];
    end
    
end
close(wb)

%%
Qc = X(:, end);
I = X(:,1);
M = X(:,2);
T = X(:,3);

plot(Qc, M, '.')
figure
plot(Qc, T, '.')


function g = ModeloLinear(xss,u)

nx = length(xss);
nu = length(u);
X = sym('X',[1 nx]).';
U = sym('U',[1 nu]).';

dXdt = derivadas(0, X, U);

Apos = jacobian(dXdt,X);
Bpos = jacobian(dXdt,U);

Apos = double(subs(Apos, [X; U], [xss'; u']));
Bpos = double(subs(Bpos, [X; U], [xss'; u']));

s = tf('s');
g = (s*eye(size(Apos,1)) - Apos)\(Bpos);

%C = eye(nx);
%statespace = ss(Apos, Bpos, C, 0);
end
