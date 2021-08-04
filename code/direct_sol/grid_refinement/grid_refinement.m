% grid_refinement.m

clear  variables; 
close  all;
clc;

%% Memory paths
addpath('../../utils');
addpath('.././config');
addpath('.././results');
addpath('../');

%% Parameters for grid refinement
tol = 1e-6;
h_min = 1;

%% Loading the corresponding configuration file
test_1_3;
rk1_prob = prob;
rk1_ode_method = ode_method;
rk1_options = options;
rk1_results_name = results_name;

test_1_4;
rk2_prob = prob;
rk2_ode_method = ode_method;
rk2_options = options;
rk2_results_name = results_name;

%% Grid refinement loop
err = tol+1;
while norm(err) > tol
    % Solving the control problems with fmincon
    tic;
    rk1_prob_sol = fmincon(@rk1_prob.F_sol,rk1_prob.z_0,[],[],[],[],rk1_prob.lb,rk1_prob.ub,@rk1_prob.nonlcon,rk1_options);
    rk1_duration_time = toc;
    fprintf('Duration time for solving the Problem: %4.2f [min]\n',rk1_duration_time/60);
    
    % Calculating the local errors on the same grid
    err = zeros(rk1_prob.N-1,1);
    for i = 1:rk1_prob.N-1
        h = rk1_prob.t(i+1)-rk1_prob.t(i);
        eta = rk1_prob_sol(i,1:4) + h * rk2_ode_method(@rk2_prob.f,[rk1_prob.t(i);rk1_prob.t(i+1)],rk1_prob_sol(i,:),2,4);
        err(i,1) = norm(rk1_prob_sol(i,1:4) - eta);
    end
    
    % Checking whether tollerance has been achieved
    if max(err) <= tol
        return; % Reduction of the maximum error
    end
%     if max(err)/min(err) <= tol
%         return; % Reduction of the uniform distribution of local errors
%     end
    
    % Calculate new grid and grid parameters
    rk1_prob.z_0 = rk1_prob_sol;
    i=1;
    j=1;
    N = rk1_prob.N;
    while j < N-2
        h = rk1_prob.t(i+1)-rk1_prob.t(i);
        if (err(j) > tol) && (norm(h) > h_min)
            rk1_prob.z_0 = [rk1_prob.z_0(1:i,:);rk1_prob.z_0(i,:);rk1_prob.z_0(i+1:end,:)];
            rk1_prob.lb = [rk1_prob.lb(1:i,:);rk1_prob.lb(i,:);rk1_prob.lb(i+1:end,:)];
            rk1_prob.ub = [rk1_prob.ub(1:i,:);rk1_prob.ub(i,:);rk1_prob.ub(i+1:end,:)];
            rk1_prob.t = [rk1_prob.t(1:i),(rk1_prob.t(i)+rk1_prob.t(i+1))/2,rk1_prob.t(i+1:end)];
            rk2_prob.t = [rk2_prob.t(1:i),(rk2_prob.t(i)+rk2_prob.t(i+1))/2,rk2_prob.t(i+1:end)];
            rk1_prob.N = rk1_prob.N + 1;
            rk2_prob.N = rk2_prob.N + 1;
            i = i + 2;
            j = j + 1;
        else
            i = i + 1;
            j = j + 1;
        end
    end
end

% Save the new grid
fprintf('Saving the grid ...\n');
writematrix(rk1_prob.t,strcat('./results/',rk1_results_name,'_grid.txt'));
writematrix(rk2_prob.t,strcat('./results/',rk2_results_name,'_grid.txt'));
fprintf('All done!\n');