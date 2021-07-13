classdef MaximalRangeFlight
    % Using a class for this provides a cleaner interface
    properties
        %% Optimalsteuerungsproblem für einen Airbus A380-800
        t_0 = 0;            % Anfangszeitpunkt in [s]
        t_f = 1800;         % Endzeitpunkt in [s]

        X_0 = [   0;         % h_0 in [m]
              27;         % gamma_0 in [rad]  
                 0;         % x_0 in [m]
              100];         % v_0 in [m/s]

        X_T = [10668;        % h_t in [m]
                  0];       % gamma_t  in [rad]

        alpha = 1.247015;   % Parameter zur Berechung der Luftdichte
        beta = 0.000104;
        g = 9.81;           % Erdbeschleunigung in [N/kg]
        C_D_0 = 0.032;      % Nullluftwiderstandsbeiwert (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
        e = 0.8;            % Oswaldfaktor (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
        F = 845;            % Wirksame Fläche / von der Luft angeströmte Fläche in [m^2]
        AR = 7.5;           % Flügelstreckung in [] (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
        k, % = 1/(pi*e*AR);    % Faktor für Berechnung des Luftwiderstandsbeiwertes
        m = 276800;         % Gewicht des Flugzeugs in [kg]
        q_max = 44154;      % Maximaler Staudruck in [N/m^2]
        T_min = 0;          % Minimale Schubkraft in [N]
        T_max = 1260000;    % Maximale Schubkraft in [N]
        C_L_min = 0.0;      % Minimaler Auftriebsbeiwert in []
        C_L_max = 1.48;     % Maximaler Auftriebsbeiwert in []
        %
        %% Problemgrößen
        n_x = 4;
        n_u = 2;
        n_c = 0;
        n_s = 1;
        n_psi = 6;
        %
        %% Grid
        N; % Anzahl an Diskretisierungen
        t;
        %
        %% Implicit or explicit iterative methods to approximate the solution of a ODE
        ode_method;
        %
        %%
        z_0;
        lb;
        ub;
    end
    methods
        %% Constructor: Set initial conditions of the direct solver
        function obj = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,ode_method)
            obj.N = N;
            obj.t = linspace(obj.t_0,obj.t_f,obj.N+1);
            obj.ode_method = ode_method;
            % Start vector
            obj.z_0 = zeros(obj.n_x+obj.n_u*(obj.N+1),1);
            obj.z_0(1:obj.n_x) = [h_0;gamma_0;x_0;v_0];
            for i = 0:N
                obj.z_0((i*obj.n_u)+obj.n_x+1:(i+1)*obj.n_u+obj.n_x) = [T_0;C_L_0];
            end
            % box conditions
            obj.lb = zeros(obj.n_x+obj.n_u*(obj.N+1),1);
            obj.ub = zeros(obj.n_x+obj.n_u*(obj.N+1),1);
            obj.lb(1:obj.n_x) = [0;-90;0;0];% eventuelle Beschränkungen an h,gamma,x,v
            obj.ub(1:obj.n_x) = [15000;90;Inf;350];
            for i = 0:N
                obj.lb((i*obj.n_u)+obj.n_x+1:(i+1)*obj.n_u+obj.n_x) = [obj.T_min;obj.C_L_min];
                obj.ub((i*obj.n_u)+obj.n_x+1:(i+1)*obj.n_u+obj.n_x) = [obj.T_max;obj.C_L_max];
            end
            %
            obj.k = 1/(pi*obj.e*obj.AR);
        end
        
        %% Differential equation
        function X = f(obj,tt,x,u) 
            X = [                                                                                  x(4)*sind(x(2));
                                       (1/(2*obj.m*x(4))) * (obj.F*u(2)*(x(4))^2*obj.alpha*exp(-obj.beta*x(1)) - 2*obj.m*obj.g*cosd(x(2)));
                                                                                                   x(4)*cosd(x(2));
                 (1/(2*obj.m)) * (2*u(1) + (-obj.C_D_0 - obj.k*(u(2))^2)*obj.F*(x(4))^2*obj.alpha*exp(-obj.beta*x(1)) - 2*obj.m*obj.g*sind(x(2)))];
        end
        
        % Funktionen für die Nebenbedingungen
        % Zielfunktional
        function T = F_sol(obj,z)
            x = obj.ode_method(@obj.f,obj.t,z,obj.N,obj.n_x,obj.n_u)';
            x = x(:,end);
            %
            T = -(x(3)-z(3));
        end
        % Ungleichungsnebenbedingungen
        function g_array = G(obj,z,x)
            g_array = zeros((obj.n_c+obj.n_s)*(obj.N+1),1);
            for i = 0:obj.N-1
                g_array((i*obj.n_s)+obj.n_c*(obj.N+1)+1) = 0.5 * obj.alpha * exp(-obj.beta*x(1,i+1)) * (x(4,i+1))^2 - obj.q_max;
            end
        end
        % Gleichungsnebenbedingungen
        function h_array = H(obj,z,x)
            x = x(:,end);
            %
            h_array = [z(1)-obj.X_0(1);
                       z(2)-obj.X_0(2);
                       z(3)-obj.X_0(3);
                       z(4)-obj.X_0(4);
                       x(1)-obj.X_T(1);
                       x(2)-obj.X_T(2)];
        end

        %% Funtion für fmincon
        % Nichtlineare Beschränkungen
        function [c,ceq] = nonlcon(obj,z)
            x = obj.ode_method(@obj.f,obj.t,z,obj.N,obj.n_x,obj.n_u)';
            c = obj.G(z,x); 
            ceq = obj.H(z,x);
        end
    end
end