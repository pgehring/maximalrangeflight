% test_1_6.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_6';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
N = 5;

direkt_sol = readmatrix(strcat('test_1_1','.txt'))';
eta_mehrsv = zeros(8,N);
for k = 0:N-1
    eta_mehrsv(:,k+1) = [direkt_sol(1:4,k+1);5.10402542449782;1145.48666826056;-1;-2.60344266831659];
end

eta_mehrsv = [-0.000461697366532953,0.986114562416963,3.94764693225948,7.71668739410971,11.0625721620951;-0.0846487572569101,0.271832230846204,0.869521184399733,0.842654452877248,0.976943535547555;4.40587449745183e-17,307.290839190894,620.813036849306,919.679653888999,1212.15142003352;99.9978160790220,104.680947420473,101.064854423120,98.5519710076505,96.4818863792967;-2.65200000524707,-2.74740160969817,-2.79123668783286,-2.79226931489557,-2.79323154274863;-1383.47152361742,-891.400583555400,-15.8162930540916,959.031785341558,1924.40581959055;-0.999999999999994,-0.999999999999996,-0.999999999999997,-0.999999999999998,-1.00000000000000;-18.9352714810131,-4.61674473792164,3.74850152514068,5.65595239889774,4.45341656466738];

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]
X_T = direkt_sol(1:2,N);
       
params = [       0,... % t_0:   Anfangszeitpunkt in [s]
               N*3,... % t_f:   Endzeitpunkt in [s]
          1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            276800,... % m:     Leergewicht des A380 in [kg]
             44154];   % q_max: Maximaler Staudruck in [N/m^2]  

prob = MaximalRangeFlightIndirect(eta_mehrsv,X_0,X_T,params);

%% Boxbeschränkungen
h_min         = 1e-12;
AbsTol        =  1e-8;
RelTol        =  1e-8;

% h_min         = 1e-8;
% AbsTol        =  1e-3;
% RelTol        =  1e-3;

StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 100;
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

% ode_method = @ode45;
ode_method = @ode23s;
shooting_methods = shooting_methods(ode_method,h_min,lb,ub,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag);
shooting_method = @shooting_methods.Mehrfachschiessverfahren;