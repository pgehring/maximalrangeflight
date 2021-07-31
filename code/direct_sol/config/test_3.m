
% mit t_0 = 1400 N = 400
h_0 = 20;     % in [m]
gamma_0 = 9; % in [Grad]  
x_0 = 6000;    % in [m]
v_0 = 90;      % in [m/s]
T_0 = 1259999;  % in [N]
C_L_0 = 1.47;   % in []

N = 400;%40;         % Anzahl an Diskretisierungen
ode_methods = ode_methods();

%Objekt der Problemklasse erhalten
% prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.irk);
prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.explicit_euler);

%Versuch1

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',4000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-8,'StepTolerance',1e-14);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-9,'StepTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',5000e+03,'MaxIterations',4.0e+05);

options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e+03,'MaxIterations',4.0e+05);

% options = optimset('Display','iter','MaxFunEvals' ,10000000,'MaxIter',4.0e+05); % Innere Punkte Verfahren
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');