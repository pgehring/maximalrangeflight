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
        n_psi = 8;
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
            obj.t = linspace(obj.t_0,obj.t_f,obj.N);
            obj.ode_method = ode_method;
            % Start vector
            obj.z_0 = zeros(obj.N,obj.n_x+obj.n_u);
            obj.z_0(1,:) = [h_0,gamma_0,x_0,v_0,T_0,C_L_0];
            % box conditions
            obj.lb = zeros(obj.N,obj.n_x+obj.n_u);
            obj.ub = zeros(obj.N,obj.n_x+obj.n_u);
            for i = 1:N
%                 obj.lb(i,:) = [0,-90,0,0,obj.T_min,obj.C_L_min]; % eventuelle Beschränkungen an h,gamma,x,v
%                 obj.ub(i,:) = [Inf,90,Inf,Inf,obj.T_max,obj.C_L_max];
                obj.lb(i,:) = [-Inf,-Inf,-Inf,-Inf,obj.T_min,obj.C_L_min]; % eventuelle Beschränkungen an h,gamma,x,v
                obj.ub(i,:) = [Inf,Inf,Inf,Inf,obj.T_max,obj.C_L_max];
            end
            %
            obj.k = 1/(pi*obj.e*obj.AR);
        end
        
        %% Differential equation
        function X = f(obj,t,z)
            X = [                                                                                                                 z(4)*sind(z(2)),...
                                              (1/(2*obj.m*z(4))) * (obj.F*z(6)*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*cosd(z(2))),...
                                                                                                                                  z(4)*cosd(z(2)),...
                 (1/(2*obj.m)) * (2*z(5) + (-obj.C_D_0 - obj.k*(z(6))^2)*obj.F*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*sind(z(2)))];
        end
        
        % Funktionen für die Nebenbedingungen
        % Zielfunktional
        function T = F_sol(obj,z)
            T = -(z(end,3)-z(1,3));
        end
        % Ungleichungsnebenbedingungen
        function g_array = G(obj,z)
            g_array = zeros(1,obj.n_s*(obj.N+1));
            for i = 0:obj.N-1
                g_array((i*obj.n_s)+1) = 0.5 * obj.alpha * exp(-obj.beta*z(i+1,1)) * (z(i+1,4))^2 - obj.q_max;
            end  
        end
        % Gleichungsnebenbedingungen
        function h_array = H(obj,z)
            h_array = zeros(1,obj.n_x*obj.N+obj.n_psi);
            x = obj.ode_method(@obj.f,obj.t,z,obj.N,obj.n_x);
            for i = 0:obj.N-2
                h_array((obj.n_x*i+1):(obj.n_x*i+obj.n_x)) = z(i+1,1:obj.n_x) + x(i+1,:) - z(i+2,1:obj.n_x);
            end
            
            
            %
            h_array((obj.n_x*obj.N+1):(obj.n_x*obj.N+obj.n_psi)) = [  z(1,1)-obj.X_0(1),...
                                                                      z(1,2)-obj.X_0(2),...
                                                                      z(1,3)-obj.X_0(3),...
                                                                      z(1,4)-obj.X_0(4),...
                                                                    z(end,1)-obj.X_T(1),...
                                                                    z(end,2)-obj.X_T(2),...
                                                                                      0,...
                                                                                      0];
        end

        %% Funtion für fmincon
        % Nichtlineare Beschränkungen
        function [c,ceq] = nonlcon(obj,z)
            c = obj.G(z); 
            ceq = obj.H(z);
        end
    end
end