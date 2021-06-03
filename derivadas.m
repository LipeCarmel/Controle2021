function [F] = derivadas(t,X,u)
%X = [I M T Tc D0 D1]

%Qi = 108;
%Qc = 471.6;
Qi = u(1);
Qc = u(2);

Qs = 459;
Qm = 378;
Qt = Qi + Qs + Qm;

Ad = 2.142*10^17;
Ed = 14897;

Ap = 3.816*10^10;
Ep = 3557;

At = 4.50*10^12;
Et = 843;

fi = 0.6;
Mm = 104.14;
deltHr = 16700;
hA = 2.52*10^5;

pCp = 360;
pcCpc = 966.3;


V = 3000;
Vc = 3312.4;

If = 0.5888;
Mf = 8.6981;

Tf = 330;
Tcf = 295;


kd = Ad*exp(-Ed/X(3));
kp = Ap*exp(-Ep/X(3));
kt = At*exp(-Et/X(3));

P = (2*fi*kd*X(1)/kt)^0.5;

F(1,1) = (Qi*If - Qt*X(1))/V - kd*X(1);
F(2,1) = (Qm*Mf - Qt*X(2))/V - kp*X(2)*P;
F(3,1) = Qt*(Tf - X(3))/V + deltHr/pCp*kp*X(2)*P - hA*(X(3)-X(4))/(pCp*V);
F(4,1) = Qc*(Tcf - X(4))/V + hA*(X(3) - X(4))/(pcCpc*Vc);
F(5,1) = 0.5*kt*P^2 - Qt*X(5)/V;
F(6,1) = Mm*kp*X(2)*P - Qt*X(6)/V;
