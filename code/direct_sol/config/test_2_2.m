% Test 2.2

%% Testparameter
N = 100;       % Anzahl an Diskretisierungen

h_sp = 9000;    % in [m]
gamma_sp = 5;   % in [Grad]  
x_sp = 800000;  % in [m]
v_sp = 250;     % in [m/s]
T_sp = 1259999; % in [N]
C_L_sp = 1.4;   % in []

t_f = 1200;    % Neue Endzeit t_f
h_0 = 100;     % Neue Starthoehe

%% Lösungsmethode der ODE und Objekt der Problemklasse erhalten
ode_methods = ode_methods();
prob = MaximalRangeFlight(h_sp,gamma_sp,x_sp,v_sp,T_sp,C_L_sp,N,@ode_methods.explicit_euler);

%% Optionen für fmincon von Matlab
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',2000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-9,'StepTolerance',1e-15);

%% Speicher Parameter
results_name = 'test_2_2';

%% Falls Daten geladen werden möchten
% test = readmatrix(strcat('./results/',results_name,'.txt'));