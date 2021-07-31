% Test 1.2

%% Testparameter
N = 100 % Anzahl an Diskretisierungen

h_0 = 9000;     % in [m]
gamma_0 = 5; % in [Grad]  
x_0 = 800000;    % in [m]
v_0 = 250;      % in [m/s]
T_0 = 1259999;  % in [N]
C_L_0 = 1.4;   % in []

t_f = 1200; % Neue Endzeit t_f

%% Lösungsmethode der ODE und Objekt der Problemklasse erhalten
ode_methods = ode_methods();
prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.explicit_euler);

%% Optionen für fmincon von Matlab
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',2000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-9,'StepTolerance',1e-15);

%% Speicher Parameter
results_name = 'test_1_1';

%% Falls Daten geladen werden möchten
% test = readmatrix(strcat('./results/',results_name,'.txt'));