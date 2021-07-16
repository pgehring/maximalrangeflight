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

addpath('../utils')
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
t = linspace(prob.t_0,prob.t_f,prob.N);
titles = [  "Flughoehe", "Anstellwinkel", ...
            "Zurueckgelegte Streckte", "Geschwindigkeit" , ...
            "Steuerung 1: Schub", "Steuerung 2: Auftriebsbeiwert"];
labels = [  "$h_{sol}\,in\,[m]$", "$\gamma_{sol}\,in\,[^{\circ}]$", ...
            "$x_{sol}\,in\,[m]$", "$v_{sol}\,in\,[\frac{m}{s}]$", ...
            "$T_{sol}\,in\,[N]$", "$C_{L_{sol}}\,in\,[1]$"];
plotter = Plotter();
plotter.plot_fmincon(t, prob_solution, titles, labels, [3, 1, 5, 2, 4, 6])

% plotter.axes(5).LineWidth = 2;
% plotter.axes(6).LineWidth = 2;