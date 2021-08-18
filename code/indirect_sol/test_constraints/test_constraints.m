% test_constraints.m
% Description:
%   Test if a state reaches the maximum back pressure condition of the
%   direct solution.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Goetz, Felix

clear  variables; 
close  all;
clc;

%% Load data
test_results_name = 'test_1_1';
run(strcat('../../direct_sol/config/',test_results_name));
data = readmatrix(strcat('../../direct_sol/results/',test_results_name,'.txt'));

%% Test if a state reaches the maximum back pressure condition
qmax_q_diff = zeros(N,3);
nof_exc = 0;
for i = 1:N
    q = (prob.alpha*exp(-prob.beta*data(i,1))*data(i,4)^2)/2;
    qmax_q_diff(i,:) = [prob.q_max,q,prob.q_max-q];
    if (q > prob.q_max)
        nof_exc = nof_exc + 1;
    end
end
fprintf('Number of times the maximum back pressure condition is exceeded: %3.0f \n',nof_exc);