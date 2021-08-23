% test_1_7.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_7';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
N = 3;
direkt_sol = readmatrix(strcat('test_1_1','.txt'))';


z_0 = [1.42507652388271,0.475382556030107,314.165065162526,107.494903804113,7.80602331112732,-1457.36355077900,-1,1.18612398450025e-05];
z_0 = [1.427991284507836;0.464716901115848;3.031184444261582e+02;1.095857118235467e+02;64.016866089723800;-1.527037711974739e+03;-1;-1.322771359854147e+02];
z_0 = [1.427991284532843;0.464716901138619;3.031184198633471e+02;1.095857107648142e+02;64.007933687835290;-1.527039943022338e+03;-1;-1.322673444184185e+02];

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)
X_0 = direkt_sol(1:4,2);
X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]
X_T = direkt_sol(1:2,3);

params = [       0,... % t_0:   Anfangszeitpunkt in [s]
                 3,... % t_f:   Endzeitpunkt in [s]
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
h_min         = 1e-12;
AbsTol        =  1e-8;
RelTol        =  1e-8;
StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 25;
% flag = 'SensDGL';
flag = 'FinitDiff';

% ode_method = @ode45;
ode_method = @ode23s;
shooting_methods = shooting_methods(ode_method,h_min,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag);
shooting_method = @shooting_methods.Einfachschiessverfahren;