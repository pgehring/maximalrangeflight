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
addpath('../direct_sol/results/');

%% Loading the corresponding configuration file
configs = ["test_1_3"];
solutions = {};
for i = 1:length(configs)
    % Load configuration file
    run(configs(i))
    % Solving the control problem with a shooting method
    try
        tic;
        [t_sol,prob_sol,eta,i_sol,Norm_F] = shooting_method( prob.tspan,...
                                                               prob.z_0,...
                                                                @prob.G,...
                                                              @prob.G_Z,...
                                                                @prob.R,...
                                                            @prob.R_Z_0,...
                                                            @prob.R_Z_f);
        duration_time = toc;
        %
        fprintf('Duration time for solving the Problem: %4.2f [min]\n',duration_time/60);
        fprintf('Required iterations: %4.2f \n',i_sol);
        if isfile(strcat('./results/',results_name,'_fmincon.txt'))
            fmincon_sol = readmatrix(strcat('./results/',results_name,'_fmincon.txt'));
            fmincon_sol(end+1,:) = [duration_time,i_sol,Norm_F(end)];
            writematrix(fmincon_sol,strcat('./results/',results_name,'_fmincon.txt'));
        else
            fmincon_sol = [duration_time,i_sol,Norm_F(end)];
            writematrix(fmincon_sol,strcat('./results/',results_name,'_fmincon.txt'));
        end
        %
        solutions{i} = {results_name,prob,prob_sol,shooting_methods.options};
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

labels_lambda = ["$\lambda_1$ in $1$","$\lambda_2$ in $1$",...
                 "$\lambda_3$ in $1$","$\lambda_4$ in $1$"];
frame_prop_lambda = [0.5,0.5,0.5,0.5];
line_style_lambda = ["b-","b-","b-","b-"];
order_lambda = [1,2,3,4];      

for j=1:length(solutions)
    solution = solutions{j};
    [results_name, prob, prob_sol, options] = solution{:};
    writematrix(prob_sol(1,:),strcat('./results/',results_name,'_vec.txt'));
    
    % Plot solution
    prob_sol_steu = prob.sol_func(prob_sol);
    fig = plotter.plot_fmincon(t_sol,prob_sol_steu,results_name,titles,labels,order,frame_prop,line_style);
    fig_lambda = plotter.plot_lambda(t_sol,prob_sol(:,5:8),results_name,labels_lambda,order_lambda,frame_prop_lambda,line_style_lambda);

    % Save the graphics
    fprintf('Saving the graphics ...\n');
    print(fig,'-dpdf','-r600',strcat('./results/',results_name,'.pdf'));
    savefig(fig,strcat('./results/',results_name,'.fig'));
    print(fig_lambda,'-dpdf','-r600',strcat('./results/',results_name,'_lambda.pdf'));
    savefig(fig_lambda,strcat('./results/',results_name,'_lambda.fig'));
    
    figure(3)
    plot(linspace(1,i_sol,length(Norm_F)),Norm_F);
    xlabel('Iterationen')
    ylabel('$\Vert F \Vert$')

    % Save the data
    fprintf('Saving the data ...\n');
    writematrix(prob_sol,strcat('./results/',results_name,'.txt'));
    writematrix(prob_sol_steu,strcat('./results/',results_name,'_steu.txt'));
end
fprintf('All done!\n');