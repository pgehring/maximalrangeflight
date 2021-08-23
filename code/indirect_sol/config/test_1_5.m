% test_1_5.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_5';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
z_0 = [ 10,...         % h_start in [m]
       0.27,...        % gamma_start in [Grad]  
         0,...         % x_start in [m]
        100,...        % v_start in [m/s]
          -25,...        % lambda_1
          400,...        % lambda_2
          10,...        % lambda_3
          -850]';          % lambda_4   
% z_0 = readmatrix(strcat('./results/',results_name,'_vec.txt'));

% z_0 = [ 5000,...         % h_start in [m]
%        25,...        % gamma_start in [Grad]  
%          10000,...         % x_start in [m]
%          300,...        % v_start in [m/s]
%           1000,...        % lambda_1
%           800000,...        % lambda_2
%           33000,...        % lambda_3
%           -20000]';          % lambda_4   
% z_0 = readmatrix(strcat('./results/',results_name,'_vec.txt'));
% 
% z_0 = [ 0,...         % h_start in [m]
%         0.27,...        % gamma_start in [Grad]  
%          0,...         % x_start in [m]
%          100,...        % v_start in [m/s]
%           -0.01,...        % lambda_1
%           -0.01,...        % lambda_2
%           -0.001,...        % lambda_3
%           -0.1]';          % lambda_4 
      
N = 5;
direkt_sol = readmatrix(strcat('test_1_1','.txt'))';
N_direkt_sol = size(direkt_sol,2);
diff = floor(N_direkt_sol/(N+1));
eta_mehrsv = zeros(length(z_0),N);
for k = 0:N-1
    eta_mehrsv(:,k+1) = [direkt_sol(1:4,(k*diff)+1);0.1;0.2;0.3;-0.1];
end

for k = 0:N-1
    eta_mehrsv(:,k+1) = [0;0.27;0;100;7.842025122992616;1.145486668260562e+03;-1;-6.242151406637181];
end

eta_mehrsv = [0,0.400000000000000,0.800000000000000,1.20000000000000,1.42799137556770;0.270000000000000,0.300000000000000,0.350000000000000,0.400000000000000,0.464716984022552;0,100,200,300,303.026938398665;100,102,104,106,109.581877190856;7.84202512299262,7.84202512299262,7.84202512299262,7.84202512299262,7.84202512299262;1145.48666826056,1145.48666826056,1145.48666826056,1145.48666826056,1145.48666826056;-1,-1,-1,-1,-1;-6.24215140663718,-6.24215140663718,-6.24215140663718,-6.24215140663718,-6.24215140663718];

eta_mehrsv = [3.71906766523057e-05,0.281930125003813,0.538704506094773,0.773604410063626,1.01738014232200;0.306709300554805,0.283836574682685,0.260709832868446,0.237112176761945,0.307092690757542;3.89630017686922e-18,60.6989507242118,122.791721126324,186.269436454574,251.005084690610;100.000157635733,102.330128915403,104.645715735755,106.946204013129,108.449509829136;5.10409756404554,5.10419953858747,5.10430193110777,5.10439586675646,5.10359107052715;1145.48666826056,819.690986161739,487.531118790968,150.210617056310,-190.485891260659;-1.00000000000000,-1.00000000000000,-1.00000000000000,-1.00000000000000,-1.00000000000000;-2.60353790153616,-2.60337110738145,-2.39527050837494,-1.99223391173045,-1.34048330604936];

eta_mehrsv = [-2.10158282196810e-05,0.281866381137973,0.538669523032691,0.773612012239903,1.01739881448357;0.297346521895342,0.268744621159702,0.243623303488853,0.221853641930065,0.297537177262935;-7.35466959318824e-18,60.6989561444505,122.791737015852,186.269469121966,251.005167112717;100.000171788594,102.330156665986,104.645756286674,106.946256523043,108.449651462932;5.10403595418053,5.10413792563131,5.10424031513867,5.10433424800055,5.10352959411690;1145.48666826056,819.695069330560,487.539278967791,150.222868743851,-190.469657240963;-0.999999999999997,-0.999999999999998,-0.999999999999999,-0.999999999999999,-1;-2.60345651029442,-2.60329264520415,-2.39519882502211,-1.99217172098314,-1.34044241607144];

eta_mehrsv = [-3.09224272039186e-05,0.281855513953035,0.538663552497432,0.773613355800593,1.01740221951113;0.295752922731381,0.266175886620184,0.240715083187094,0.219256558996444,0.295910781453914;4.23581212610560e-18,60.6989572415536,122.791740354926,186.269475918577,251.005182852220;100.000174205437,102.330161406759,104.645763213340,106.946265484923,108.449675533296;5.10402542449782,5.10412739542293,5.10422978441795,5.10432371680594,5.10351908714837;1145.48666826056,819.695765900468,487.540669546824,150.224955041300,-190.466893030637;-1.00000000000000,-1.00000000000000,-1.00000000000000,-1.00000000000000,-1;-2.60344266831659,-2.60327930053318,-2.39518663191069,-1.99216114171554,-1.34043546536608];

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
maxit = 2;
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