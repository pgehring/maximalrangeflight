% MaximalRangeFlight.m:

% Klasse des Optimalsteuerungsproblem für einen Airbus A380-800

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

classdef MaximalRangeFlight
    properties
        %% Parameter für das Optimalsteuerungsproblem
        t_0 = 0;            % Anfangszeitpunkt in [s]
        t_f = 1400;         % Endzeitpunkt in [s]

        X_0 = [   0;        % h_0 in [m]
               0.27;        % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
                  0;        % x_0 in [m]
                100];       % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

        X_T = [10668;       % h_t in [m]
                   0];      % gamma_t  in [Grad]

        alpha = 1.247015;   % Parameter zur Berechung der Luftdichte in []
        beta = 0.000104;
        g = 9.81;           % Erdbeschleunigung in [N/s^2]
        C_D_0 = 0.032;      % Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
        e = 0.8;            % Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
        F = 845;            % Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
        AR = 7.5;           % Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
        k;                  % Faktor für Berechnung des Luftwiderstandsbeiwertes in []
        m = 500000;%570000;% Startgewicht normal 276800;         % Leergewicht des A380 in [kg]
        q_max = 44154;      % Maximaler Staudruck in [N/m^2]
        T_min = 0;          % Minimale Schubkraft in [N]
        T_max = 1260000;    % Maximale Schubkraft in [N]
        C_L_min = 0.0;      % Minimaler Auftriebsbeiwert in []
        C_L_max = 1.48;     % Maximaler Auftriebsbeiwert in []

        %% Problemgrößen
        n_x = 4;            % Größe des Zustandsvektors
        n_u = 2;            % Größe des Steuervektors
        n_c = 0;            % Größe des gemischten Steuer- und Zustandsbeschränkungsvektors
        n_s = 1;            % Größe des reinen Zustandsbeschränkungsvektors
        n_psi = 8;          % Grüße des Randbedingungsvektors

        %% Variablen für die Berechnung
        N;                  % Variable für die anzahl an Diskretisierungen
        t;                  % Variable für den Zeitenvektor
        ode_method;         % Implicit or explicit iterative methods to approximate the solution of a ODE
        z_0;                % Array für die Startwerte für fmincon
        lb;                 % lower bounds for fmincon
        ub;                 % upper bounds for fmincon
    end
    methods
        %% Constructor: Set initial conditions of the direct solver
        function obj = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,ode_method, z_0)
            % grid size
            obj.N = N;
            obj.t = linspace(obj.t_0,obj.t_f,obj.N);
            
            obj.ode_method = ode_method;
            
            % start vector
            if nargin > 8
                obj.z_0 = z_0;
            else
                obj.z_0 = [h_0,gamma_0,x_0,v_0,T_0,C_L_0].*ones(obj.N,obj.n_x+obj.n_u);
            end
            
            % box conditions
            obj.lb = zeros(obj.N,obj.n_x+obj.n_u);
            obj.ub = zeros(obj.N,obj.n_x+obj.n_u);
            for i = 1:N
%                 obj.lb(i,:) = [          0,...
%                                        -90,...
%                                          0,...
%                                          0,...
%                                  obj.T_min,...
%                                obj.C_L_min];
%                 obj.ub(i,:) = [      13100,... % Maximale Flughöhe A380 in [m]
%                                         90,...
%                                   15200000,... % Maximale Reichweite A380 in [m]
%                                        266,... % Maximale Geschwindigkeit A380 in [m/s]
%                                  obj.T_max,...
%                                obj.C_L_max];
                obj.lb(i,:) = [          -inf,...
                                       -inf,...
                                         -inf,...
                                         -inf,...
                                 obj.T_min,...
                               obj.C_L_min];
                obj.ub(i,:) = [      inf,... % Maximale Flughöhe A380 in [m]
                                        inf,...
                                  inf,... % Maximale Reichweite A380 in [m]
                                       inf,... % Maximale Geschwindigkeit A380 in [m/s]
                                 obj.T_max,...
                               obj.C_L_max];
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
        
        %% Funktionen für die Nebenbedingungen
        % Zielfunktional
        function T = F_sol(obj,z)
            T = -(z(end,3)-obj.X_0(3));
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
            
%             % Test mit ode von Matlab
%             y = z(1,1:obj.n_x);
%             for i = 0:obj.N-2
%                 func =@(t,x) obj.f(t,[x',z(i+1,obj.n_x:end)])';
%                 dt = abs((obj.t(i+1)-obj.t(i+2))/10);
%                 tspan = [obj.t(i+1):dt:obj.t(i+2)];
%                 sol = ode23s(func,tspan,y);
%                 y = sol.y(:,end)';
%                 h_array((obj.n_x*i+1):(obj.n_x*i+obj.n_x)) = y - z(i+2,1:obj.n_x);
%             end
            
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