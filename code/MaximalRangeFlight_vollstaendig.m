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

%% Optimalsteuerungsproblem für einen Airbus A380-800
% Parameter:
t_0 = 0;            % Anfangszeitpunkt in [s]
t_f = 1800;         % Endzeitpunkt in [s]

X0 = [   0;         % h_0 in [m]
      0.27;         % gamma_0 in [rad]  
         0;         % x_0 in [m]
      100];         % v_0 in [m/s]
         
XT = [10668;        % h_t in [m]
          0];       % gamma_t  in [rad]
      
alpha = 1.247015;   % Parameter zur Berechung der Luftdichte
beta = 0.000104;
g = 9.81;           % Erdbeschleunigung in [N/kg]
C_D_0 = 0.032;      % Nullluftwiderstandsbeiwert (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
e = 0.8;            % Oswaldfaktor (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
F = 845;            % Wirksame Fläche / von der Luft angeströmte Fläche in [m^2]
AR = 7.5;           % Flügelstreckung in [] (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
k = 1/(pi*e*AR);    % Faktor für Berechnung des Luftwiderstandsbeiwertes
m = 276800;         % Gewicht des Flugzeugs in [kg]
q_max = 44154;      % Maximaler Staudruck in [N/m^2]
T_min = 0;          % Minimale Schubkraft in [N]
T_max = 1260000;    % Maximale Schubkraft in [N]
C_L_min = 0.0;      % Minimaler Auftriebsbeiwert in []
C_L_max = 1.48;     % Maximaler Auftriebsbeiwert in []










%% Problembeschreibung
% Problemgrößen
n_x = 4;
n_u = 2;
n_c = 0;
n_s = 1;
n_psi = 6;
% Gitter
N = 30; % Anzahl an Diskretisierungen
t = linspace(t_0,t_f,N+1);
% Anfangs- und Endbedingungen
x_start = X0; % Anfangsgeschwindigkeit
x_end = XT;

%% Optimierungsproblem
% Startvektor
x0 = [   0;         % h_0 in [m]
      -45;         % gamma_0 in [rad]  
         0;         % x_0 in [m]
      100];         % v_0 in [m/s]
u0 = [T_max;        % T
      C_L_min];       % C_L
z0 = zeros((n_x+n_u)*(N+1),1);
for i = 0:N
    z0((i*n_x)+1:(i+1)*n_x) = x0;
    z0((i*n_u)+n_x*(N+1)+1:(i+1)*n_u+n_x*(N+1)) = u0;
end
% Boxbeschränkungen
lb = zeros((n_x+n_u)*(N+1),1);
ub = zeros((n_x+n_u)*(N+1),1);
for i = 0:N
    lb((i*n_x)+1:(i+1)*n_x) = [0;-90;0;0]; % eventuelle Beschränkungen an h,gamma,x,v
    lb((i*n_u)+n_x*(N+1)+1:(i+1)*n_u+n_x*(N+1)) = [T_min;C_L_min];
    ub((i*n_x)+1:(i+1)*n_x) = [15000;90;1e8;350];
    ub((i*n_u)+n_x*(N+1)+1:(i+1)*n_u+n_x*(N+1)) = [T_max;C_L_max];
end
% Zielfunktional minimieren
F_sol =@(z) -(z(n_x*N+3)-z(3));
% Nichtlineare Beschränkungen
nonlcon =@(z)my_nonlcon(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max);
%
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',60.000000e+03);
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%
Sol = fmincon(F_sol,z0,[],[],[],[],lb,ub,nonlcon,options);









%% Plot der Lösungen
for i = 1:N+1
    h_sol(i) = Sol(((i-1)*n_x)+1);
    gamma_sol(i) = Sol(((i-1)*n_x)+2);
    x_sol(i) = Sol(((i-1)*n_x)+3);
    v_sol(i) = Sol(((i-1)*n_x)+4);
    T_sol(i) = Sol(n_x*(N+1)+((i-1)*n_u)+1);
    C_L_sol(i) = Sol(n_x*(N+1)+((i-1)*n_u)+2);
end
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








%% Funktionen für die Nebenbedingungen
% Ungleichungsnebenbedingungen
function g_array = Ungleichungsnebenbedingungen(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max)
g_array = zeros((n_c+n_s)*(N+1),1);

h = @(i)      z(n_x*i+1);
gamma = @(i)  z(n_x*i+2);
x = @(i)      z(n_x*i+3);
v = @(i)      z(n_x*i+4);
T = @(i)      z((i*n_u)+n_x*(N+1)+1);
C_L = @(i)    z((i*n_u)+n_x*(N+1)+2);
              
for i = 0:N-1
    g_array((i*n_s)+n_c*(N+1)+1) = 0.5 * alpha * exp(-beta*h(i)) * (v(i))^2 - q_max;
end  
end
% Gleichungsnebenbedingungen
function h_array = Gleichungsnebenbedingungen(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max)
h_array = zeros(n_x*N+n_psi,1);

h = @(i)      z(n_x*i+1);
gamma = @(i)  z(n_x*i+2);
x = @(i)      z(n_x*i+3);
v = @(i)      z(n_x*i+4);
T = @(i)      z((i*n_u)+n_x*(N+1)+1);
C_L = @(i)    z((i*n_u)+n_x*(N+1)+2);

% DGL
f =@(t,i) [                                                                                   v(i)*sind(gamma(i));
                                (1/(2*m*v(i))) * (F*C_L(i)*(v(i))^2*alpha*exp(-beta*h(i)) - 2*m*g*cosd(gamma(i)));
                                                                                              v(i)*cosd(gamma(i));
           (1/(2*m)) * (2*T(i) + (-C_D_0 - k*(C_L(i))^2)*F*(v(i))^2*alpha*exp(-beta*h(i)) - 2*m*g*sind(gamma(i)))];

for i = 0:N-1
    h_array((i*n_x)+1:(i+1)*n_x) = z((i*n_x)+1:(i+1)*n_x) + (t(i+2)-t(i+1))*f(t(i+1),i) - z(((i+1)*n_x)+1:(i+2)*n_x);
end  

h_array(n_x*N+1:n_x*N+n_psi) = [    z(1)-x_start(1);
                                    z(2)-x_start(2);
                                    z(3)-x_start(3);
                                    z(4)-x_start(4);
                                z(n_x*N+1)-x_end(1);
                                z(n_x*N+2)-x_end(2)];
end

% Funtion für fmincon
function [c,ceq] = my_nonlcon(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max)
c = Ungleichungsnebenbedingungen(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max); 
ceq = Gleichungsnebenbedingungen(z,t,N,x_start,x_end,n_x,n_u,n_c,n_s,n_psi,alpha,beta,g,C_D_0,e,F,AR,k,m,q_max,T_min,T_max,C_L_min,C_L_max);
end
