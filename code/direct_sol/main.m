% main.m:
% Task:
%   Model of a 2-dimensional flight of an aircraft in the x-h plane, where 
%   the lift coefficient C_L and the thrust T can be controlled. A maximum 
%   dynamic pressure q must not be exceeded. The objective is to control 
%   the flight of an aircraft from a given initial position in such a way 
%   that a given cruise altitude of 10668 m is reached, the angle of attack 
%   there is 0 degrees and the range is maximum.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

clear  variables
close  all
clc

%% Memory paths
addpath('../utils')
addpath('./config')
addpath('./results')

%% Loading the corresponding configuration file
test_1_1
% test_1_2

% test_2_1
% test_2_1

% test_3_1
% test_3_2

% test_4_1
% test_4_2

% test_5_1

% test_6_1

%% Solving the control problem with fmincon
tic;
prob_sol = fmincon(@prob.F_sol,prob.z_0,[],[],[],[],prob.lb,prob.ub,@prob.nonlcon,options);
duration_time = toc;
fprintf('Duration time for solving the Problem: %4.2f [min]\n', duration_time/60);

%% Plot the solution and saving the results
t = linspace(prob.t_0,prob.t_f,prob.N);
titles = ["Flughoehe","Anstellwinkel","Zurueckgelegte Streckte",...
          "Geschwindigkeit","Steuerung 1: Schub",...
          "Steuerung 2: Auftriebsbeiwert"];
labels = ["$h_{sol}$ in $[m]$","$\gamma_{sol}$ in $[^{\circ}]$",...
          "$x_{sol}$ in $[m]$","$v_{sol}$ in $[\frac{m}{s}]$",...
          "$T_{sol}$ in $[N]$","$C_{L_{sol}}$ in $[1]$"];
frame_prop = [0.5,0.5,0.5,0.5,2,2];
line_style = ["b-","b-","b-","b-","r-","r-"];
order = [3,1,5,2,4,6];
plotter = Plotter();
fig = plotter.plot_fmincon(t,prob_sol,titles,labels,order,frame_prop,line_style);

% Save the graphics
fprintf('Saving the graphics ...\n');
savefig(fig,strcat('./results/',results_name,'.fig'))
saveas(fig,strcat('./results/',results_name,'.png'))
saveas(fig,strcat('./results/',results_name,'.svg'))

% Save the data
fprintf('Saving the data ...\n');
writematrix(prob_sol,strcat('./results/',results_name,'.txt'))

fprintf('All done!\n');