% test_1_1.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_3';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
% z_0 = [ 10,...         % h_start in [m]
%        0.27,...        % gamma_start in [Grad]  
%          0,...         % x_start in [m]
%         100,...        % v_start in [m/s]
%           -25,...        % lambda_1
%           400,...        % lambda_2
%           10,...        % lambda_3
%           -850];          % lambda_4   
% z_0 = readmatrix(strcat('./results/',results_name,'_vec.txt'));

z_0 = [ 5000,...         % h_start in [m]
       25,...        % gamma_start in [Grad]  
         10000,...         % x_start in [m]
         300,...        % v_start in [m/s]
          1000,...        % lambda_1
          800000,...        % lambda_2
          33000,...        % lambda_3
          -20000];          % lambda_4   
% z_0 = readmatrix(strcat('./results/',results_name,'_vec.txt'));

z_0 = [ 0,...         % h_start in [m]
        0.27,...        % gamma_start in [Grad]  
         0,...         % x_start in [m]
         100,...        % v_start in [m/s]
          -0.01,...        % lambda_1
          -0.001,...        % lambda_2
          -0.1,...        % lambda_3
          -0.01];          % lambda_4   
      
% 1.) Versuch
% z_0 = [7.500213289569822;0.269999999999992;0;100;-30.294341505117618;3.054837129407352e+02;7.250234618526340;-8.456721736033198e+02]
% [1.335411545171215e-21;0.270000000000000;0;100;-29.488316939766570;3.889531824497242e+02;-1;-8.519023137701304e+02]

N = 3;
direkt_sol = readmatrix(strcat('test_1_1','.txt'))';
z_0 = [direkt_sol(1:4,1);1000;8000000;33000;-2000]; %;1;1;1;1];   

% z_0(5:8,1) = [11986768.9505666,761594731.477802,6.67046151566364e-06,-13021165.6489014]';

X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]
X_T = direkt_sol(1:2,2);

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

prob = maximal_range_flight(z_0,X_0,X_T,params);

%% Boxbeschränkungen
h_min         = 1e-12;
AbsTol        =  1e-8;
RelTol        =  1e-8;

StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 1;
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

ode_method = @ode45;
ode_method = @ode15s;
shooting_methods = shooting_methods(ode_method,h_min,lb,ub,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag);
shooting_method = @shooting_methods.Einfachschiessverfahren;