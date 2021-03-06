classdef ReatorPolimer < handle
    properties
        % Valor default
        Mw = 0;
        
        % Par?metros:
        % Vaz?es
        Qm = 378;
        % Par?metros cineticos
        Ad = 2.142*10^17;
        Ed = 14897;
        
        Ap = 3.816*10^10;
        Ep = 3557;
        
        At = 4.50*10^12;
        Et = 843;
        
        deltHr = 16700;
        hA = 2.52*10^5;
        
        % Propriedades de pol?meros
        fi = 0.6;
        Mm = 104.14;
        
        % Propiedades do fluido refrigerante
        pCp = 360;
        pcCpc = 966.3;
        
        % Volumes
        V = 3000;
        Vc = 3312.4;
        
        % Concentra??es de alimenta??o
        If = 0.5888;
        Mf = 8.6981;
        
        % Temperaturas de alimenta??o
        Tf = 330;
        Tcf = 295;
    end
    methods
        function obj = ReatorPolimer()
        end
        function [F] = derivadas(obj, t, X, u)
            %X = [I M T Tc D0 D1]
            %Qi = 108;
            %Qc = 471.6;
            Qi = u(1);
            Qc = u(2);
            
            Qs = 1.5*obj.Qm - Qi;
            Qt = Qi + Qs + obj.Qm;
            
            obj.Mw = X(6)/X(5);
            
            kd = obj.Ad*exp(-obj.Ed/X(3));
            kp = obj.Ap*exp(-obj.Ep/X(3));
            kt = obj.At*exp(-obj.Et/X(3));
            
            P = (2*obj.fi*kd*X(1)/kt)^0.5;
            
            F(1,1) = (Qi*obj.If - Qt*X(1))/obj.V - kd*X(1);
            F(2,1) = (obj.Qm*obj.Mf - Qt*X(2))/obj.V - kp*X(2)*P;
            F(3,1) = Qt*(obj.Tf - X(3))/obj.V + obj.deltHr/obj.pCp*kp*X(2)*P - obj.hA*(X(3)-X(4))/(obj.pCp*obj.V);
            F(4,1) = Qc*(obj.Tcf - X(4))/obj.V + obj.hA*(X(3) - X(4))/(obj.pcCpc*obj.Vc);
            F(5,1) = 0.5*kt*P^2 - Qt*X(5)/obj.V;
            F(6,1) = obj.Mm*kp*X(2)*P - Qt*X(6)/obj.V;
        end
        function n = viscosidade(obj)
           n = 0.0012*obj.Mw^0.71; 
        end
    end
    methods(Static)
        function visc = vetor_viscosidade(D0, D1)
            visc = 0.0012*(D1./D0).^0.71;
        end
    end
end