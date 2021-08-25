% main.m:
% Description:
%   Model of a 2-dimensional flight of an aircraft in the x-h plane, where 
%   the lift coefficient C_L and the thrust T can be controlled. A maximum 
%   dynamic pressure q must not be exceeded. The objective is to control 
%   the flight of an aircraft from a given initial position in such a way 
%   that a given cruise altitude of 10668 m is reached, the angle of attack 
%   there is 0 degrees and the range is maximum.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Goetz, Felix

clear  variables; 
close  all;
clc;

%% Setup workspace
addpath('../utils');
addpath('./config');
addpath('./results');

%% Loading the corresponding configuration file
configs = ["test_0_1","test_0_2"];
% configs = ["test_2_1","test_2_2"];
% configs = ["test_1_1"];
% configs = ["test_3_1","test_3_2"];
% configs = ["test_4_1","test_4_2"];
solutions = {};
for i = 1:length(configs)
    % Load configuration file
    run(configs(i))
    % Solving the control problem with fmincon
    try
        tic;
        [prob_sol,fval,exitflag,output] = fmincon(@prob.F_sol,prob.z_0,[],[],[],[],prob.lb,prob.ub,@prob.nonlcon,options);
        duration_time = toc;
        %
        fprintf('Duration time for solving the Problem: %4.2f [min]\n',duration_time/60);
        fprintf('Required iterations: %4.2f \n',output.iterations);
        fprintf('Required function evaluations: %4.2f \n',output.funcCount);
        if isfile(strcat('./results/',results_name,'_fmincon.txt'))
            fmincon_sol = readmatrix(strcat('./results/',results_name,'_fmincon.txt'));
            fmincon_sol(end+1,:) = [fval,exitflag,duration_time,output.iterations,output.funcCount,output.firstorderopt];
            writematrix(fmincon_sol,strcat('./results/',results_name,'_fmincon.txt'));
        else
            fmincon_sol = [fval,exitflag,duration_time,output.iterations,output.funcCount,output.firstorderopt];
            writematrix(fmincon_sol,strcat('./results/',results_name,'_fmincon.txt'));
        end
        %
        solutions{i} = {results_name, prob, prob_sol, options};
    catch ME
        fprintf('Error occured while solving control problem with config %s\nContinuing with next file..\n', configs(i))
        disp(ME.message)
        continue % with next config
    end
end

%% Plot the solution and saving the results
plotter = Plotter();
titles = [{'Flugh\"ohe'},{'Anstellwinkel'},{'Zur\"uckgelegte Streckte'},...
          {'Geschwindigkeit'},{'\textbf{Steuerung 1: Schub}'},...
          {'\textbf{Steuerung 2: Auftriebsbeiwert}'}];
labels = ["$h_{sol}$ in $m$","$\gamma_{sol}$ in $^{\circ}$",...
          "$x_{sol}$ in $m$","$v_{sol}$ in $\frac{m}{s}$",...
          "$T_{sol}$ in $N$","$C_{L_{sol}}$ in $1$"];
frame_prop = [0.5,0.5,0.5,0.5,2,2];
line_style = ["b-","b-","b-","b-","r-","r-"];
order = [3,1,5,2,4,6];

for j=1:length(solutions)
    solution = solutions{j};
    [results_name, prob, prob_sol, options] = solution{:};
    % Plot solution
    fig = plotter.plot_fmincon(prob.t,prob_sol,results_name,titles,labels,order,frame_prop,line_style);
    plotter.plot_staudruck(prob.t,prob_sol,prob,results_name);
    % Save the graphics
    fprintf('Saving the graphics ...\n');
    print(fig,'-dpdf','-r600',strcat('./results/',results_name,'.pdf'));
    savefig(fig,strcat('./results/',results_name,'.fig'));
    %
    fmincon_gesamt = readmatrix(strcat('./results/',results_name,'_fmincon.txt'));
    Aufwand(1,1) = strcat("Funktionswert: ",string(fmincon_gesamt(end,1)));
    Aufwand(2,1) = strcat("Exitflag: ",string(fmincon_gesamt(end,2)));
    Aufwand(3,1) = strcat("Berechnungsdauer in Minuten: ",string(sum(fmincon_gesamt(:,3))/60));
    Aufwand(4,1) = strcat("Anzahl Funktionsauswertungen: ",string(sum(fmincon_gesamt(:,5))));
    Aufwand(5,1) = strcat("Anzahl Iterationen: ",string(sum(fmincon_gesamt(:,4))));
    Aufwand(6,1) = strcat("Optimalitaetsmass: ",string(fmincon_gesamt(end,6)));
    writematrix(Aufwand,strcat('./results/',results_name,'_Aufwand.txt'));
    % Save the data
    fprintf('Saving the data ...\n');
    writematrix(prob_sol,strcat('./results/',results_name,'.txt'));
end
fprintf('All done!\n');