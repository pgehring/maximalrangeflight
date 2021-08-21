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

MODEL = "d" % model switch  d: direct_model
            %               i: indirect_model
            

%% solution parameters
N = 100;

t0 = 0;
tf = 600;
t = linspace(t0, tf, N);


switch MODEL
    case "i"
        X0 = [  0;           % h_0 in [m]
                0.27;           % gamma_0 in [rad]
                0;           % x_0 in [m]
                100;        % v_0 in [m/s] 
                -1;        
                -1;           
                -1;           
                1];          
%         X0(2) = X0(2),*180/pi;
        
        odemethods = {@ode45, @ode23s, @euler_expl};
        ode_labels = ["ode45", "ode23s", "euler expl"];
        ode_styles = ["--r", ":m", "--b", "c"];
        ref_method_idx = 3;
        
        func = @indirect_model;
        names = ["methods_plot_i_h","methods_plot_i_gamma", "methods_plot_i_x", "methods_plot_i_v", "methods_plot_i_l1", "methods_plot_i_l2", "methods_plot_i_l3", "methods_plot_i_l4"];
        titles = [{'Flugh\"ohe'},{'Anstellwinkel'},{'Zur\"ueckgelegte Streckte'},{'Geschwindigkeit'},{''},{''},{''},{''}];
        labels = ["$h_{sol}$ in $[m]$","$\gamma_{sol}$ in $[^{\circ}]$","$x_{sol}$ in $[m]$","$v_{sol}$ in $[\frac{m}{s}]$", "$\lambda_1$", "$\lambda_2$", "$\lambda_3$", "$\lambda_4$"];
    case "d"
        func = @direct_model;
        names = ["methods_plot_d_h","methods_plot_d_gamma", "methods_plot_d_x", "methods_plot_d_v"];
        titles = [{'Flugh\"ohe'},{'Anstellwinkel'},{'Zur\"ueckgelegte Streckte'},{'Geschwindigkeit'}];
        labels = ["$h_{sol}$ in $[m]$","$\gamma_{sol}$ in $[^{\circ}]$","$x_{sol}$ in $[m]$","$v_{sol}$ in $[\frac{m}{s}]$"];
              
        odemethods = {@ode45, @ode23s, @euler_expl, @euler_impl, @irk};
        ode_labels = ["ode45", "ode23s", "euler expl", "euler impl", "radau-2a"];
        ode_styles = ["--r", ":m", "--b", "-.k", "c"];
        ref_method_idx = 3;
        
        X0 = [  0;           % h_0 in [m]
                0.27;           % gamma_0 in [rad]
                0;           % x_0 in [m]
                100];          % v_0 in [m/s] 
%         X0(2) = X0(2)*180/pi;
end

%% solve model
solutions = {};
for i=1:length(odemethods)
    solver = odemethods{i};
    opts = odeset('RelTol',1e-3,'AbsTol',1e-5);
    
    tic
    [~, X] = solver(func, t, X0, opts, params_cell);
    time_elapsed = toc;
    fprintf("solution took %fs to complete with %s\n", time_elapsed, ode_labels(i))
    solutions{i} = {X, time_elapsed, ode_labels(i)};
end

% compute maximum state deviation and time performance
ref_vector = solutions{ref_method_idx}{1};
ref_time = solutions{ref_method_idx}{2};
for i=1:length(odemethods)
    solution_vec = solutions{i}{1};
    
    % state deviation
    deviation = solution_vec - ref_vector;
    y_max = max(abs(deviation))
    
    % time performance
    time = solutions{i}{2};
    time_delta = time - ref_time
    
    solutions{i} = [solutions{i}, {y_max, time_delta}];
    
end

%% plotting
plot_position = [50, 50, 500, 350];

set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLineMarkerSize',4);
set(0, 'DefaultLineLineWidth', 0.8);
FigW=7.5;
FigH=4;
loc=["northwest","northeast","northwest","southeast","southwest","northwest","northwest","northwest"];

figures = [];
for jj=1:length(titles)
    f = figure(jj);
    f.Position = plot_position;
    figures = [figures, f ];
    set(f,'defaulttextinterpreter','latex',...
            'PaperUnits','centimeters','PaperSize',[FigW FigH],...
                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
                'Position',[0,0,FigW,FigH]);

    set(f, 'PaperPositionMode', 'manual');
    hold on
    for j = 1:length(solutions)
        disp(ode_labels(j))
        [X_plot, ~, ~, ~] = solutions{j}{:};
        
        plot(gca(), t, X_plot(:, jj), ode_styles(j))
        ylabel(sprintf("%s\n%s", string(titles(jj)), labels(jj)),'FontSize',7)
        xlabel("Zeit $t$ in $s$",'FontSize',7)
    end
    hold off
    set(findall(gcf,'-property','FontSize'),'FontSize',8);
    [lgd, objh] = legend(ode_labels, 'Interpreter','latex','Location',loc(jj),'FontSize',5,'Orientation','vertical');
    set(f, 'PaperPositionMode', 'auto');
    ax = findobj(f, 'type', 'axes');
    ax.Units = 'centimeters';
    test=ax.Position;
    ax.Position=[test(1)+0.2,test(2),test(3),test(4)];
    print(f,'-dpdf','-r600',strcat('./results/', names(jj),'.pdf'));
end

%% Workspace cleanup
diary off