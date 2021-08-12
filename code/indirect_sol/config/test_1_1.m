% test_1_1.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_1';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
z_0 = [  0,...         % h_start in [m]
       0.27,...        % gamma_start in [Grad]  
         0,...         % x_start in [m]
        100,...        % v_start in [m/s]
          10,...        % lambda_1
          2,...        % lambda_2
          1.5,...        % lambda_3
          3];          % lambda_4   

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]

params = [       0,... % t_0:   Anfangszeitpunkt in [s]
               550,... % t_f:   Endzeitpunkt in [s]
          1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            530000,... % m:     Leergewicht des A380 in [kg]
             44154];   % q_max: Maximaler Staudruck in [N/m^2]  

prob = maximal_range_flight(z_0,X_0,X_T,params);

%% Boxbeschränkungen
h_min         = 1e-10;
AbsTol        =  1e-8;
RelTol        =  1e-8;
StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 1000;
flag = 'FinitDiff';
% flag = 'SensDGL';
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

shooting_methods = shooting_methods(h_min,lb,ub,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag);
shooting_method = @shooting_methods.Einfachschiessverfahren;