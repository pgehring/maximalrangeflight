% test_1_4.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_4';

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
z_0 = [direkt_sol(1:4,1);0.1;0.2;0.3;-0.1];   

% z_0 = [0,0.270000000000017,0,100,426.263251767933,140621.198044161,-1,-151.848994695743]';

%z_0 = [-1.88577859267535e-15,0.270000000000017,0,100,455.277037798810,82311.7221922134,-1,-90.1464713728160]';

% z_0(5:8,1) = [11986768.9505666,761594731.477802,6.67046151566364e-06,-13021165.6489014]';


z_0 = [0;0.27;0;100;20;1000;35;-100];  
z_0 = [-1.42744220487496e-14;0.270000000000002;0;100;9.28257800133037;1177.90345832235;-1.00000000004258;-8.31631681339405];
z_0 = [-1.427442202482772e-14;0.270000000000002;0;100;9.282578012148559;1.177903458317930e+03;-1.000000000042580;-8.316316830859838];
z_0 = [-6.648161524994853e-15;0.270000000000001;0;100;8.032927230956808;1.145364979391497e+03;-1.000000000014971;-6.546693988736394];
z_0 = [0;0.27;0;100;8.033;1.145364979391497e+03;-1.0;-6.546693988736394];
z_0 = [-1.682949869594878e-25;0.270000000000000;0;100;8.033000002653100;1.145364979409142e+03;-1;-6.546693992751067];
z_0 = [0;0.27;0;100;8.034;1.145364979409142e+03;-1;-6.546693992751067];
z_0 = [2.827983158496073e-26;0.270000000000000;0;100;8.034000000089295;1.145364979437482e+03;-1;-6.546693992802615];
z_0 = [2.827983158496073e-26;0.270000000000000;0;100;8.035;1.145364979437482e+03;-1;-6.546693992802615];
z_0 = [-6.646763760486976e-26;0.270000000000000;0;100;8.035000000066768;1.145364979437813e+03;-1;-6.546693992903565];
z_0 = [-6.646763760486976e-26;0.270000000000000;0;100;8.0355000000066768;1.145364979437813e+03;-1;-6.546693992903565];
z_0 = [3.259482265368129e-17;0.270000000000000;0;100;7.959301222095917;1.145324348345823e+03;-1;-6.418165113009415];
z_0 = [-5.946684281312355e-17;0.270000000000000;0;100;7.845449920619672;1.145484908321891e+03;-1;-6.247080117320514];
z_0 = [-5.940095237168155e-17;0.270000000000000;0;100;7.845095507427610;1.145483661170902e+03;-1;-6.246624032085068];
z_0 = [-5.856728396249574e-17;0.270000000000000;0;100;7.843401267164286;1.145487560553867e+03;-1;-6.244771286919990];
z_0 = [-5.856728397125160e-17;0.270000000000000;0;100;7.843401266990457;1.145487560554003e+03;-1;-6.244771286811942];
z_0 = [-5.732440736775094e-17;0.270000000000000;0;100;7.842025122992616;1.145486668260562e+03;-1;-6.242151406637181];
z_0 = [-5.732440736775094e-17;0.270000000000000;0;100;7.842025122992616;1.145486668260562e+03;-1;-6.242151406637181];




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
maxit = 50;
% flag = 'SensDGL';
flag = 'FinitDiff';
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
shooting_method = @shooting_methods.Einfachschiessverfahren;