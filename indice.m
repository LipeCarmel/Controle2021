classdef indice
    methods(Static)
        function indice = ISE(e, Ts)
            f = e.^2;
            indice = Ts*trapz(f);
        end
        function indice = IAE(e, Ts)
            f = abs(e);
            indice = Ts*trapz(f);
        end
        function indice = ITSE(e, Ts)
            t = 0 : Ts : Ts*(length(e) - 1);
            f = t.*e.^2;
            indice = Ts*trapz(f);
        end
        function indice = ITAE(e, Ts)
            t = 0 : Ts : Ts*(length(e) - 1);
            f = abs(t.*e);
            indice = Ts*trapz(f);
        end
        % Teste
        function resultado = teste_de_validacao()
            n = 6;
            e = ones(1, n);
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