% irk.m:
% implicit Runge-Kutta procedure for a butcher Tableau
% Shema: 
%   [ c | A;
%       | b]
% Input:    f       - function  handle  for rhs of the DGL
%           tspan   - boundaries of the interval [t_start , t_end]
%           A       - A from Butcher Tableau
%           b       - b from Butcher Tableau
%           c       - c from Butcher Tableau
%           y0      - initial  values y_0
%           N       - number of equidistant time N N
% Output:   t       - discretized  time of tspan
%           y       - solution of the ODE
% Date:         18.05.2021
% Author:       Felix Goetz

function y = irk(f,t,z,N,n_x,n_u)

% RADAU-2A-Verfahren
A = [5/12, -1/12;
      3/4,   1/4];
b = [3/4; 1/4];
c = [1/3; 1];

% s-stufiges Runge-Kutta Verfahren
s = length(c); 
% K-Werte von Runge-Kutta
K = zeros(length(z(1:n_x)),s);
% Lösungsmatrix
y = zeros(N+1,length(z(1:n_x)));
y(1,:) = z(1:n_x);
% Parameter für fsolve
options = optimoptions('fsolve','Display','none','OptimalityTolerance',1e-8);
% Berechenen der Werte für jeden Zeitpunkt
for n = 1:N
    % Schrittweite
    h = (t(n+1)-t(n));
    % Steuerung u_n
    u(1) = z((n*n_u)+n_x+1);
    u(2) = z((n*n_u)+n_x+1);
    % 
    K_New =@(K) stufenableitung(f,y(n,:),s,t(n),h,A,c,K,u);
    K = fsolve(K_New,K,options);
    % Runge-Kutta Schritt
    y(n+1,:) = y(n,:) + h*b'*K';
end
end

function K_new = stufenableitung(f,y,s,t,h,A,c,K,u)
K_new = K;
for i = 1:s
    K_new(:,i) = K(:,i) - f(t+c(i)*h, y+h*A(i,:)*K',u);
end
end
