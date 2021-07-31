% Test 3.1

%% Testparameter
N = 100;       % Anzahl an Diskretisierungen

h_sp = 20;      % in [m]
gamma_sp = 9;   % in [Grad]  
x_sp = 6000;    % in [m]
v_sp = 90;      % in [m/s]
T_sp = 1259999; % in [N]
C_L_sp = 1.47;  % in []

t_f = 1200;    % Neue Endzeit t_f
m = 500000;    % Neues Gewicht

%% Lösungsmethode der ODE und Objekt der Problemklasse erhalten
ode_methods = ode_methods();
prob = MaximalRangeFlight(h_sp,gamma_sp,x_sp,v_sp,T_sp,C_L_sp,N,@ode_methods.explicit_euler);

%% Optionen für fmincon von Matlab
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',2000.0e+03,'MaxIterations',4.0e+05);

%% Speicher Parameter
results_name = 'test_3_1';

%% Falls Daten geladen werden möchten
% test = readmatrix(strcat('./results/',results_name,'.txt'));