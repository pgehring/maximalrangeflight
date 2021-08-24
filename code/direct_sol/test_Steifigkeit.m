clear  variables; 
close  all;
clc;

%% Setup workspace
addpath('../utils');
addpath('./config');
addpath('./results');

%%
results_name = 'test_1_1';

run(results_name)
prob.z_0 = readmatrix(strcat('./results/',results_name,'.txt'));

%[0,86.7561354086353,0,0.497330169101140;-2.92972839148699e-05,0.0487880895822648,0,0.00366812421791290;0,-49.7330169167982,0,0.867561353969753;0.000293377318110555,-8.51077688244328,0,-0.0564187150136779]

eig_val = zeros(N,1);
nof_exc = 0;
for i = 1:N
    Jac = G_Z(prob.z_0(i,:),prob);
    Jac_app = CalcJacobian(@prob.f,prob.z_0(i,:));
    Eigenvalues = real(eig(Jac));
    
    abs_eig = abs(Eigenvalues(2:4));
    testt(i,1) = max(abs_eig)/min(abs_eig);
    
end
x = linspace(0,prob.t_f-1,N);
plot(x,testt(:,1),'r');

function jacobian = CalcJacobian(fun,x)
    x = x(:);
    x(2) = x(2)*pi/180;
    t = fun(1,x);
    z = t(1:4);
    n=numel(x);
    m=numel(z);
    jacobian = zeros(m,m);
    h= 1e-6;
    for i=1:m
        % Einheitsvektor
        unit_vector = zeros(n,1);
        unit_vector(i) = 1;

        % Neues x
        x_tmp = x+h*unit_vector;

        % Berechne die Jacobi-matrix
        t = fun(1,x_tmp);
        jacobian(:,i) = (t(1:4)-z)./h;
    end
end

function Z = G_Z(z,prob)
    z(2) = z(2)*pi/180;
    % Synthese-Steuerung
    T = z(5);
    C_L = z(6);
    %
    J_1_2 = + z(4)*cos(z(2));
    J_1_4 = + sin(z(2));
    %
    J_2_1 = - (prob.F*prob.alpha*prob.beta*exp(-prob.beta*z(1))*z(4)*C_L)/(2*prob.m);
    J_2_2 = + (prob.g*sin(z(2)))/z(4);
    J_2_4 = + (prob.F*prob.alpha*exp(-prob.beta*z(1))*C_L)/(2*prob.m)...
            + (prob.g*cos(z(2)))/(z(4)^2);
    %
    J_3_2 = - z(4)*sin(z(2));
    J_3_4 = + cos(z(2));
    %
    J_4_1 = + ((prob.C_D_0+prob.k*C_L^2)*prob.F*prob.alpha*prob.beta*exp(-prob.beta*z(1))*z(4)^2)/(2*prob.m);
    J_4_2 = - prob.g*cos(z(2));
    J_4_4 = - ((prob.C_D_0+prob.k*C_L^2)*prob.F*prob.alpha*exp(-prob.beta*z(1))*z(4))/(prob.m);
    %
    Z = [    0,J_1_2,    0,J_1_4;
         J_2_1,J_2_2,    0,J_2_4;
             0,J_3_2,    0,J_3_4;
         J_4_1,J_4_2,    0,J_4_4];
end