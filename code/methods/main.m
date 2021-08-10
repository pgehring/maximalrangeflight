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

clear; 
close  all;
clc;

%% Setup workspace
addpath('../utils');

diary off
if exist('console.log', 'file')
    delete('console.log')
end

%% model parameters
params = [1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            500000,... % m:     Leergewicht des A380 in [kg]
             44154,... % q_max: Maximaler Staudruck in [N/m^2]  
           1259999,... % T_start in [N]
              1.47];   % C_L_start in [1]
params_cell = num2cell(params);

%% solution parameters
N = 100;

t0 = 0;
tf = 1200;
t = linspace(t0, tf, N);

X0 = [  0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] 
X0(2) = X0(2)*180/pi;

odemethods = {@ode45, @ode23s, @euler_expl, @euler_impl, @irk};
ode_labels = ["ode45", "ode23s", "euler expl", "euler impl", "radau-2a"];
ode_styles = ["--r", ":m", "--b", "-.k", "c"];
%% solve model
solutions = {};
for i=1:length(odemethods)
    solver = odemethods{i};
    opts = odeset('RelTol',1e-2,'AbsTol',1e-5);
    
    tic
    [~, X] = solver(@model, t, X0, [], params_cell);
    toc
%     fprintf("time elapsed")
    solutions{i} = {X};
end

%% plotting
titles = ["Flughoehe","Anstellwinkel","Zurueckgelegte Streckte",...
          "Geschwindigkeit"];
labels = ["$h_{sol}$ in $[m]$","$\gamma_{sol}$ in $[^{\circ}]$",...
          "$x_{sol}$ in $[m]$","$v_{sol}$ in $[\frac{m}{s}]$"];
plot_position = [50, 50, 800, 400];

set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
      
figures = [];
for jj=1:length(titles)
    f = figure(jj);
    f.Position = plot_position;
    figures = [figures, f ];
    
    hold on
    for j = 1:length(solutions)
        disp(ode_labels(jj))
        [X_plot] = solutions{j}{:};
        
        plot(gca(), t, X_plot(:, jj), ode_styles(j))
        ylabel(sprintf("%s\n%s", titles(jj), labels(jj)))
        xlabel("Zeit $t$ in $s$")
    end
    hold off
    legend(ode_labels, 'Interpreter','latex')
end

%% Workspace cleanup
diary off