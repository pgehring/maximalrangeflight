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
h_0 = 2131;     % in [m]
gamma_0 = 0.27; % in [Grad]  
x_0 = 12312;    % in [m]
v_0 = 100;      % in [m/s]
T_0 = 1260000;  % in [N]
C_L_0 = 1.48;   % in []

N = 35;         % Anzahl an Diskretisierungen
ode_methods = ode_methods();

%% Objekt der Problemklasse erhalten
% prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.irk);
prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.explicit_euler);

%% Fmincon
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',300.0e+03,'MaxIterations',4.0e+05);
% options = optimset('Display','iter','MaxFunEvals' ,190000,'MaxIter',4.0e+05); % Innere Punkte Verfahren
%
tic;
prob_solution = fmincon(@prob.F_sol,prob.z_0,[],[],[],[],prob.lb,prob.ub,@prob.nonlcon,options);
duration_time = toc

%% Plot der Lösungen
fig = figure(1);
fig.WindowState = 'maximized'; % Vollbild
t = linspace(prob.t_0,prob.t_f,prob.N);

subplot(2,3,1);
plot(t,prob_solution(:,3),'b-');
title('Zurückgelegte Streckte');
ylabel('x_{sol} in [m]');
xlabel('t in [s]');

subplot(2,3,2);
plot(t,prob_solution(:,1),'b-');
title('Höhe');
ylabel('h_{sol} in [m]');
xlabel('t in [s]');

T_plot = subplot(2,3,3);
plot(t,prob_solution(:,5),'r-');
title('Steuerung 1: Schub');
ylabel('T_{sol} in [N]');
xlabel('t in [s]');
T_plot.LineWidth = 2;

subplot(2,3,4);
plot(t,prob_solution(:,2),'b-');
title('Anstellwinkel');
ylabel('gamma_{sol} in [Grad]');
xlabel('t in [s]');

subplot(2,3,5);
plot(t,prob_solution(:,4),'b-');
title('Geschwindigkeit');
ylabel('v_{sol} in [m/s]');
xlabel('t in [s]');

C_L_plot = subplot(2,3,6);
plot(t,prob_solution(:,6),'r-');
title('Steuerung 2: Auftriebsbeiwert');
ylabel('C_{L_{sol}} in []');
xlabel('t in [s]');
C_L_plot.LineWidth = 2;