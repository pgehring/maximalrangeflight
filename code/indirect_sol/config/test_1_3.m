% test_1_3.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_3';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
z_0 = [-5.732440736775094e-17;0.270000000000000;0;100;7.842025122992616;1.145486668260562e+03;-1;-6.242151406637181];
z_0 = [0.028019454974903;0.632328123755776;-6.883525059503436e-16;1.000134561046371e+02;11.408733692221102;1.144903630842599e+03;-0.903886575198159;-7.358110154112469];
z_0 = [0.027950743602661;1.001868462525390;-5.184395014501834e-15;1.000160983759926e+02;16.619900713120284;1.144903630842599e+03;-0.904314993678545;-29.464517741834314];

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]

params = [       0,... % t_0:   Anfangszeitpunkt in [s]
               300,... % t_f:   Endzeitpunkt in [s]
          1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            276800,... % m:     Leergewicht des A380 in [kg]
             44154];   % q_max: Maximaler Staudruck in [N/m^2]  

prob = MaximalRangeFlightIndirect(z_0,X_0,X_T,params);

%% Boxbeschränkungen
h_min         = 1e-6;
AbsTol        =  1e-8;
RelTol        =  1e-8;
StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 20;
flag = 'SensDGL';
% flag = 'FinitDiff';

ode_method = @ode45;
% ode_method = @ode23s;
shooting_methods = shooting_methods(ode_method,h_min,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag);
shooting_method = @shooting_methods.Einfachschiessverfahren;