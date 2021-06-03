classdef indice
    methods(Static)
        function integral = ISE(e, Ts)
            f = e.^2;
            integral = Ts*trapz(f);
        end
        function integral = IAE(e, Ts)
            f = abs(e);
            integral = Ts*trapz(f);
        end
        function integral = ITSE(e, Ts)
            t = [0 : Ts : Ts*(length(e) - 1)]';
            f = t.*(e.^2);
            integral = Ts*trapz(f);
        end
        function integral = ITAE(e, Ts)
            t = [0 : Ts : Ts*(length(e) - 1)]';
            f = abs(t.*e);
            integral = Ts*trapz(f);
        end
        % Teste
        function resultado = teste_de_validacao()
            n = 6;
            e = ones(n, 1);
            Ts = 1;
            
            % Integral de constante
            % Resultado (n-1)*Ts, utiliza-se n-1 pois t(0) = 0
            IAE = indice.IAE(e, Ts);
            ISE = indice.ISE(e, Ts);   % Erro(k) = 1, resultado equivalente
            Ref_const = (n-1)*Ts;
            
            % Integral da reta
            % Resultado ((n-1)*Ts)^2/2, utiliza-se n-1 pois t(0) = 0
            ITAE = indice.ITAE(e, Ts);
            ITSE = indice.ITSE(e, Ts); % Erro(k) = 1, resultado equivalente
            Ref_reta = ((n-1)*Ts)^2/2;
            resultado = table(ISE, IAE, Ref_const, ITSE, ITAE, Ref_reta);
        end
    end
end