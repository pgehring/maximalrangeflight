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
% h_0 = 2131;     % in [m]
% gamma_0 = 0.27; % in [Grad]  
% x_0 = 12312;    % in [m]
% v_0 = 100;      % in [m/s]
% T_0 = 1260000;  % in [N]
% C_L_0 = 1.48;   % in []

% h_0 = 423;     % in [m]
% gamma_0 = 45; % in [Grad]  
% x_0 = 123;    % in [m]
% v_0 = 204;      % in [m/s]
% T_0 = 0;  % in [N]
% C_L_0 = 0;   % in []

% h_0 = 5000;     % in [m]
% gamma_0 = 45; % in [Grad]  
% x_0 = 50000;    % in [m]
% v_0 = 204;      % in [m/s]
% T_0 = 0;  % in [N]
% C_L_0 = 0;   % in []

% h_0 = 6000;     % in [m]
% gamma_0 = 45; % in [Grad]  
% x_0 = 60000;    % in [m]
% v_0 = 250;      % in [m/s]
% T_0 = 0;  % in [N]
% C_L_0 = 0;   % in []

% % möglicher startpunkt für active Set 
% h_0 = 0;     % in [m] 
% gamma_0 = 0.27; % in [Grad]  
% x_0 = 0;    % in [m]
% v_0 = 10;      % in [m/s]
% T_0 = 1561;  % in [N]
% C_L_0 = 0.5;   % in []

% % mit t_0 = 1200 N = 100
% h_0 = 9000;%6000;%3000;     % in [m]
% gamma_0 = 5;%45; % in [Grad]  
% x_0 = 800000;%60000;    % in [m]
% v_0 = 250;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.4;%0.1;   % in []

% % mit t_0 = 1200 N = 100 gamma_lb = -90
% h_0 = 20;     % in [m]
% gamma_0 = 9; % in [Grad]  
% x_0 = 6000;    % in [m]
% v_0 = 90;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.47;   % in []

% % mit t_0 = 1400 N=100
% h_0 = 3000;%6000;%3000;     % in [m]
% gamma_0 = -2;%45; % in [Grad]  
% x_0 = 600000;%60000;    % in [m]
% v_0 = 250;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.4;%0.1;   % in []


% % mit t_0 = 1400 N = 100 IMPLIZIT
% h_0 = 9000;%6000;%3000;     % in [m]
% gamma_0 = 5;%45; % in [Grad]  
% x_0 = 800000;%60000;    % in [m]
% v_0 = 250;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.4;%0.1;   % in []

% % mit t_0 = 1400 N = 250 gamma_lb = -90
% h_0 = 20;     % in [m]
% gamma_0 = 9; % in [Grad]  
% x_0 = 6000;    % in [m]
% v_0 = 90;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.47;   % in []

% % mit t_0 = 1400 N = 300 gamma_lb = -90
% h_0 = 20;     % in [m]
% gamma_0 = 9; % in [Grad]  
% x_0 = 6000;    % in [m]
% v_0 = 90;      % in [m/s]
% T_0 = 1259999;  % in [N]
% C_L_0 = 1.47;   % in []
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1000.0e+03);

% mit t_0 = 1400 N = 400
h_0 = 20;     % in [m]
gamma_0 = 9; % in [Grad]  
x_0 = 6000;    % in [m]
v_0 = 90;      % in [m/s]
T_0 = 1259999;  % in [N]
C_L_0 = 1.47;   % in []

N = 400;%40;         % Anzahl an Diskretisierungen
ode_methods = ode_methods();

%% Objekt der Problemklasse erhalten
% prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.irk);
prob = MaximalRangeFlight(h_0,gamma_0,x_0,v_0,T_0,C_L_0,N,@ode_methods.explicit_euler);

%Versuch1

%% Fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',4000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-8,'StepTolerance',1e-14);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-9,'StepTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',5000e+03,'MaxIterations',4.0e+05);

options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e+03,'MaxIterations',4.0e+05);

% options = optimset('Display','iter','MaxFunEvals' ,10000000,'MaxIter',4.0e+05); % Innere Punkte Verfahren
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
%
tic;
prob_solution = fmincon(@prob.F_sol,prob.z_0,[],[],[],[],prob.lb,prob.ub,@prob.nonlcon,options);
duration_time = toc

%% Plot der Lösungen
t = linspace(prob.t_0,prob.t_f,prob.N);
titles = [  "Flughoehe", "Anstellwinkel", ...
            "Zurueckgelegte Streckte", "Geschwindigkeit" , ...
            "Steuerung 1: Schub", "Steuerung 2: Auftriebsbeiwert"];
labels = [  "$h_{sol}$ in $[m]$", "$\gamma_{sol}$ in $[^{\circ}]$", ...
            "$x_{sol}$ in $[m]$", "$v_{sol}$ in $[\frac{m}{s}]$", ...
            "$T_{sol}$ in $[N]$", "$C_{L_{sol}}$ in $[1]$"];
frame_prop = [0.5,0.5,0.5,0.5,2,2];
line_style = ["b-","b-","b-","b-","r-","r-"];
plotter = Plotter();
plotter.plot_fmincon(t, prob_solution, titles, labels, [3, 1, 5, 2, 4, 6], frame_prop, line_style)

%% Automatisches Abpeichern der Daten und der Grafik


