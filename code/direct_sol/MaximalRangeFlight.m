% MaximalRangeFlight.m:

% Klasse des Optimalsteuerungsproblem f�r einen Airbus A380-800

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / G�tz, Felix

classdef MaximalRangeFlight
    properties
        %% Parameter f�r das Optimalsteuerungsproblem
        t_0 = 0;            % Anfangszeitpunkt in [s]
        t_f = 1400;         % Endzeitpunkt in [s]

        X_0 = [   0;        % h_0 in [m]
               0.27;        % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20�)
                  0;        % x_0 in [m]
                100];       % v_0 in [m/s] (Ben�tigte Startgeschwindigkeit zum Abheben)

        X_T = [10668;       % h_t in [m]
                   0];      % gamma_t  in [Grad]

        alpha = 1.247015;   % Parameter zur Berechung der Luftdichte in []
        beta = 0.000104;
        g = 9.81;           % Erdbeschleunigung in [N/s^2]
        C_D_0 = 0.032;      % Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Gr��e, Geschwindigkeit und Flugh�he in Beziehung setzt)
        e = 0.8;            % Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die �nderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Fl�gels oder Flugzeugs im Vergleich zu einem idealen Fl�gel mit demselben Seitenverh�ltnis darstellt)
        F = 845;            % Wirksame Fl�che in [m^2]                     (von der Luft angestr�mte Fl�che) 
        AR = 7.5;           % Fl�gelstreckung in []                        (aspect ratio) (Das Seitenverh�ltnis eines Fl�gels ist definiert als das Verh�ltnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragfl�cheninhalt in [m^2] = b^2 / F
        k;                  % Faktor f�r Berechnung des Luftwiderstandsbeiwertes in []
        m = 500000;%570000;% Startgewicht normal 276800;         % Leergewicht des A380 in [kg]
        q_max = 44154;      % Maximaler Staudruck in [N/m^2]
        T_min = 0;          % Minimale Schubkraft in [N]
        T_max = 1260000;    % Maximale Schubkraft in [N]
        C_L_min = 0.0;      % Minimaler Auftriebsbeiwert in []
        C_L_max = 1.48;     % Maximaler Auftriebsbeiwert in []

        %% Problemgr��en
        n_x = 4;            % Gr��e des Zustandsvektors
        n_u = 2;            % Gr��e des Steuervektors
        n_c = 0;            % Gr��e des gemischten Steuer- und Zustandsbeschr�nkungsvektors
        n_s = 1;            % Gr��e des reinen Zustandsbeschr�nkungsvektors
        n_psi = 8;          % Gr��e des Randbedingungsvektors

        %% Variablen f�r die Berechnung
        N;                  % Variable f�r die anzahl an Diskretisierungen
        t;                  % Variable f�r den Zeitenvektor
        ode_method;         % Implicit or explicit iterative methods to approximate the solution of a ODE
        z_0;                % Array f�r die Startwerte f�r fmincon
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
            obj.lb = [-inf, -inf, -inf, -inf, obj.T_min, obj.C_L_min].*ones(obj.N,obj.n_x+obj.n_u);
            obj.ub = [inf, inf, inf, inf, obj.T_max, obj.C_L_max].*ones(obj.N,obj.n_x+obj.n_u);

            obj.k = 1/(pi*obj.e*obj.AR);
        end
        
        %% Differential equation
        function X = f(obj,t,z)
            X = [                                                                                                                 z(4)*sind(z(2)),...
                                              (1/(2*obj.m*z(4))) * (obj.F*z(6)*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*cosd(z(2))),...
                                                                                                                                  z(4)*cosd(z(2)),...
                 (1/(2*obj.m)) * (2*z(5) + (-obj.C_D_0 - obj.k*(z(6))^2)*obj.F*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*sind(z(2)))];
        end
        
        % subject function
        function T = F_sol(obj,z)
            T = -(z(end,3)-obj.X_0(3));
        end
        
        %% nonlinear constraint function
        function [c,ceq] = nonlcon(obj,z)
            
            % inequality constraints G(z)
            c = 0.5 * obj.alpha * exp(-obj.beta*z(:,1)) .* (z(:,4)).^2 - obj.q_max;
            
            % equality constraints H(z)
            x = obj.ode_method(@obj.f,obj.t,z,obj.N,obj.n_x);
            
            ceq = reshape(z(obj.n_x+1:end-obj.n_x, 1:obj.n_x)', 1, []) + reshape(x(obj.n_x+1:end-obj.n_x,:)', 1, []) - reshape(z(2*obj.n_x+1:end, 1:obj.n_x)', 1, []);
            ceq((obj.n_x*obj.N+1):(obj.n_x*obj.N+obj.n_psi)) = [  z(1,1)-obj.X_0(1),...
                                                                  z(1,2)-obj.X_0(2),...
                                                                  z(1,3)-obj.X_0(3),...
                                                                  z(1,4)-obj.X_0(4),...
                                                                z(end,1)-obj.X_T(1),...
                                                                z(end,2)-obj.X_T(2),...
                                                                                  0,...
                                                                                  0];
        end
    end
end