clear all
close all
clc

u = [108, 471.6];
y0 = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];

u = [208, 471.6];
y0 = [6.6832e-2, 3.3245, 323.56, 305.17, 2.7547e-4, 16.110];


[t,y] = ode15s(@(t,y) derivadas(t,y,u), [0: .1 : 200], y0);
%%

            %X = [I M T Tc D0 D1]
figure
plot(t,y(:,1))
figure
plot(t,y(:,2))
figure
plot(t,y(:,3))
figure
plot(t,y(:,4))
figure
plot(t,(y(:,end-1)./y(:,end)).^0.71)
%%
% subplot(311)
% plot(t,y(:,end-1)./y(:,end),'k-');
% subplot(312)
% plot(t,y(:,end-1));
% subplot(313)
% plot(t,y(:,end))
