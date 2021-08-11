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

function [t,y] = irk(f,tspan,y0,opts, param)
    A = [5/12, -1/12;
      3/4,   1/4];
    b = [3/4; 1/4];
    c = [1/3; 1];
    
    % s-stufiges Runge-Kutta Verfahren
    s = length(c); 
    % K-Werte von Runge-Kutta
    K = zeros(length(y0),s);
    
    if length(tspan) == 2
        t0 = tspan(1);
        tf = tspan(2);
        h = (tf-t0)/N;
        t = t0:h:tf;
    else
        h = tspan(2)-tspan(1);
        t = tspan;
    end
    N = length(t)-1;
    y=zeros(N+1,length(y0));
    y(1,:)=y0(:)';
    
    % Parameter für fsolve
    options = optimoptions('fsolve','Display','none','Algorithm', 'levenberg-marquardt', 'OptimalityTolerance', 1e-8);
    % Berechenen der Werte für jeden Zeitpunkt
    for n = 1:N
        % 
        K_New =@(K) stufenableitung(f,y(n,:),s,t(n),h,A,c,K, param);
        K = fsolve(K_New,K,options);
        % Runge-Kutta Schritt
        y(n+1,:) = y(n,:) + h*b'*K';
    end
end

function K_new = stufenableitung(f,y,s,t,h,A,c,K, param)
    K_new = K;
    for i = 1:s
        K_new(:,i) = K(:,i) - f(t+c(i)*h, y+h*A(i,:)*K', param);
    end
end
