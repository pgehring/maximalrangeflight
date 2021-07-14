% ode_method.m:

% Klasse an verschiedenen Einschrittverfahren zur Lösung von ODE's

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

classdef ode_methods
    properties
        %% Variablen für die Berechnung
        A;
        b;
        c;
        s;
    end
    methods
        %% Constructor
        function obj = ode_methods()
            % RADAU-2A-Verfahren
            obj.A = [5/12, -1/12;
                      3/4,   1/4];
            obj.b = [3/4; 1/4];
            obj.c = [1/3; 1];
            % s-stufiges Runge-Kutta Verfahren
            obj.s = length(obj.c); 
        end
        
        %% Explicit method
        function Z = explicit_euler(obj,f,t,z,N,n_x)
            Z = zeros(N,n_x);
            % Berechenen der Werte für jeden Zeitpunkt
            for n = 1:N-1
                Z(n,:) = (t(n+1)-t(n))*f(t(n),z(n,:));
            end
        end
        
        %% Implicit methods
        function Z = irk(obj,f,t,z,N,n_x)
            % K-Werte von Runge-Kutta
            K = zeros(obj.s,n_x);
            % Lösungsmatrix
            Z = zeros(N,n_x);
            % Parameter für fsolve
            options = optimoptions('fsolve','Display','none','OptimalityTolerance',1e-8);
            % Berechenen der Werte für jeden Zeitpunkt
            for n = 1:N-1
                % Schrittweite
                h = (t(n+1)-t(n));
                % Berechnen der Vektoren k_i für i = 1,..,s
                K_New =@(K) obj.stufenableitung(f,z(n,:),t(n),h,K,n_x);
                K = fsolve(K_New,K,options);
                % Runge-Kutta Schritt
                Z(n,:) = h*obj.b'*K;
            end
        end
        
        function K_new = stufenableitung(obj,f,y,t,h,K,n_x)
            K_new = K;
            for i = 1:obj.s
                y(1:n_x) = y(1:4)+h*obj.A(1,:)*K;
                K_new(i,:) = K(i,:) - f(t+obj.c(i)*h,y);
            end
        end
        
    end
end
