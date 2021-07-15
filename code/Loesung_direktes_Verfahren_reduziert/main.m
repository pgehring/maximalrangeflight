% MaximalRangeFlight.m:

% Aufgabenstellung:
%   Modell eines 2-dimensionalen Fluges eines Flugzeugs in der x-h-Ebene,
%   bei dem der Auftriebsbeiwert und der Schub gesteuert werden kann. 
%   Dabei darf ein maximaler Staudruck nicht überschritten werden.
%   Das Ziel ist den Flug eines Flugzeuges von einer gegebenen 
%   Anfangsposition so zu steuern, dass eine vorgegebene Reisehöhe erreicht 
%   wird, der Anstellwinkel dort 0 Grad und die Reichweite maximal ist.

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

% Quellen:
%   - https://www.calculatoratoz.com/de/zero-lift-drag-coefficidet-at-minimum-required-thrust-calculator/Calc-5844
%   - https://tu-dresden.de/ing/maschinenwesen/ilr/ressourcen/dateien/tfd/studium/dateien/Flugmechanik_V.pdf?lang=de
%   - https://tu-dresden.de/ing/maschinenwesen/ilr/ressourcen/dateien/tfd/studium/dateien/Flugmechanik_U.pdf?lang=de


clear  variables
close  all
clc

%% Problembeschreibung: Optimalsteuerungsproblem für einen Airbus A380-800
h_0 = 1;         %  in [m]
gamma_0 = 30;         %  in [rad]  
x_0 = 1;         %  in [m]
v_0 = 100;         %  in [m/s]
T_0 = 600000;        % T
C_L_0 = 0.9;       % C_L
N = 10; % Anzahl an Diskretisierungen
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
plot(prob_solution)
for i = 1:N+1
    T_sol(i) = Sol(n_x+((i-1)*n_u)+1);
    C_L_sol(i) = Sol(n_x+((i-1)*n_u)+2);
end
% DGL
f =@(X,z,t,i) [                                                                                   X(4)*sind(X(2));
                                  (1/(2*m*X(4))) * (F*z((i*n_u)+n_x+2)*(X(4))^2*alpha*exp(-beta*X(1)) - 2*m*g*cosd(X(2)));
                                                                                                X(4)*cosd(X(2));
              (1/(2*m)) * (2*z((i*n_u)+n_x+1) + (-C_D_0 - k*(z((i*n_u)+n_x+2))^2)*F*(X(4))^2*alpha*exp(-beta*X(1)) - 2*m*g*sind(X(2)))];
X(:,1) = Sol(1:n_x);
for i = 0:N-1
    X(:,2+i) = X(:,1+i) + (t(i+2)-t(i+1))*f(X,Sol,t(i+1),i);
end
h_sol = X(1,:);
gamma_sol = X(2,:);
x_sol = X(3,:);
v_sol = X(4,:);

figure(1)
plot(t,x_sol,'r-',t,T_sol,'k:');
legend('x_{sol}','T_{sol}');

figure(2)
plot(t,h_sol);
legend('h_{sol}');

figure(3)
plot(t,v_sol,'b-');
legend('v_{sol}');

figure(4)
plot(t,gamma_sol,'r-',t,C_L_sol,'g*-');
legend('gamma_{sol}','C_{L_{sol}}');

% plot(t,h_sol,t,gamma_sol,t,x_sol,'r-',t,v_sol,'b-',t,T_sol,'k:',t,C_L_sol,'g*-');
% legend('h_{sol}','gamma_{sol}','x_{sol}','v_{sol}','T_{sol}','C_{L_{sol}}');







