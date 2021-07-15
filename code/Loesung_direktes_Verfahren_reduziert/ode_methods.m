classdef ode_methods
    % Using a class for this provides a cleaner interface
    properties
        A;
        b;
        c;
    end
    methods
        % Constructor
        function obj = ode_methods()
            % RADAU-2A-Verfahren
            obj.A = [5/12, -1/12;
                  3/4,   1/4];
            obj.b = [3/4; 1/4];
            obj.c = [1/3; 1];
        end

        function y = irk(obj,f,t,z,N,n_x,n_u)
            % s-stufiges Runge-Kutta Verfahren
            s = length(obj.c); 
            % K-Werte von Runge-Kutta
            K = zeros(length(z(1:n_x)),s);
            % Lösungsmatrix
            y = zeros(N+1,length(z(1:n_x)));
            y(1,:) = z(1:n_x);
            % Parameter für fsolve
            options = optimoptions('fsolve','Display','none','OptimalityTolerance',1e-8);
            % Berechenen der Werte für jeden Zeitpunkt
            for n = 1:N
                % Schrittweite
                h = (t(n+1)-t(n));
                % Steuerung u_n
                u(1) = z((n*n_u)+n_x+1);
                u(2) = z((n*n_u)+n_x+1);
                % 
                K_New =@(K) obj.stufenableitung(f,y(n,:),s,t(n),h,K,u);
                K = fsolve(K_New,K,options);
                % Runge-Kutta Schritt
                y(n+1,:) = y(n,:) + h*obj.b'*K';
            end
        end

        function K_new = stufenableitung(obj,f,y,s,t,h,K,u)
            K_new = K;
            for i = 1:s
                K_new(:,i) = K(:,i) - f(t+obj.c(i)*h, y+h*obj.A(i,:)*K',u);
            end
        end
    end
end
