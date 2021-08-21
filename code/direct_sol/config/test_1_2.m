% test_1_2.m

% Versuchsaufbau:
%   - Veränderte Endzeit: t_f = 350 [s]
%   - Explizites Euler-Verfahren für ODE
%   - SQP-Verfahren

%% Speicher Parameter
results_name = 'test_1_2';

%% Testparameter
N = 100;               % Anzahl an Diskretisierungen

z_0 = [   9000,...     % h_start in [m]
             5,...     % gamma_start in [Grad]  
        800000,...     % x_start in [m]
           250,...     % v_start in [m/s]
       1259999,...     % T_start in [N]
           1.4];       % C_L_start in []
z_0 = readmatrix(strcat('./results/',results_name,'.txt')); % Falls Daten geladen werden möchten     

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [Grad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]

params = [       0,... % t_0:   Anfangszeitpunkt in [s]
               350,... %900,... % t_f:   Endzeitpunkt in [s]
          1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            276800,... % m:     Leergewicht des A380 in [kg]
             44154];   % q_max: Maximaler Staudruck in [N/m^2]  

t = linspace(params(1),params(2),N);

%% Boxbeschränkungen
lb = [   -inf,...
         -inf,...
         -inf,...
         -inf,...
          0.0,...      % T_min:   Minimale Schubkraft in [N]
          0.0];        % C_L_min: Minimaler Auftriebsbeiwert in []
      
ub = [    inf,...
          inf,...
          inf,...
          inf,...
      1260000,...      % T_max:   Maximale Schubkraft in [N]
         1.48];        % C_L_max: Maximaler Auftriebsbeiwert in []

%% Lösungsmethode der ODE und Objekt der Problemklasse erhalten
ode_methods = ode_methods();
ode_method = @ode_methods.explicit_euler;
prob = MaximalRangeFlight(N,t,z_0,X_0,X_T,params,lb,ub,ode_method);

%% Optionen für fmincon von Matlab
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1.0e+03,'MaxIterations',4.0e+05);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',6000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-8,'UseParallel',true);