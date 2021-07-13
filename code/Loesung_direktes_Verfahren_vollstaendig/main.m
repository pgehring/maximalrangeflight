% main.m:

% Aufgabenstellung:
%   Modell eines 2-dimensionalen Fluges eines Flugzeugs in der x-h-Ebene,
%   bei dem der Auftriebsbeiwert und der Schub gesteuert werden kann. 
%   Dabei darf ein maximaler Staudruck nicht überschritten werden.
%   Das Ziel ist den Flug eines Flugzeuges von einer gegebenen 
%   Anfangsposition so zu steuern, dass eine vorgegebene Reisehöhe erreicht 
%   wird, der Anstellwinkel dort 0 Grad und die Reichweite maximal ist.

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

clear  variables
close  all
clc

%% Problembeschreibung: Optimalsteuerungsproblem für einen Airbus A380-800
h_0 = 100;         %  in [m]
gamma_0 = 0;         %  in [Grad]  
x_0 = -36;         %  in [m]
v_0 = -3;         %  in [m/s]
T_0 = 1260000;        % T
C_L_0 = 1.48;       % C_L
N = 100; % Anzahl an Diskretisierungen
ode_methods = ode_methods();

%%
prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.irk);

%% Fmincon
%
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',60.000000e+03);
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%
prob_solution = fmincon(@prob.F_sol,prob.z_0,[],[],[],[],prob.lb,prob.ub,@prob.nonlcon,options);

%% Plot der Lösungen
n_x = prob.n_x;
n_u = prob.n_u;
t = linspace(prob.t_0,prob.t_f,prob.N);
for i = 1:N+1
    h_sol(i) = prob_solution(((i-1)*n_x)+1);
    gamma_sol(i) = prob_solution(((i-1)*n_x)+2);
    x_sol(i) = prob_solution(((i-1)*n_x)+3);
    v_sol(i) = prob_solution(((i-1)*n_x)+4);
    T_sol(i) = prob_solution(n_x*(N+1)+((i-1)*n_u)+1);
    C_L_sol(i) = prob_solution(n_x*(N+1)+((i-1)*n_u)+2);
end
figure(1)
plot(t,x_sol,'r-',t,T_sol,'k:');
legend('x_{prob_solution}','T_{prob_solution}');

figure(2)
plot(t,h_sol);
legend('h_{prob_solution}');

figure(3)
plot(t,v_sol,'b-');
legend('v_{prob_solution}');

figure(4)
plot(t,gamma_sol,'r-',t,C_L_sol,'g*-');
legend('gamma_{prob_solution}','C_{L_{prob_solution}}');