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
C_D_0 = 0.032;      % Nullluftwiderstandsbeiwert (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
e = 0.8;            % Oswaldfaktor (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
F = 845;            % Wirksame Fläche / von der Luft angeströmte Fläche in [m^2]
AR = 7.5;           % Flügelstreckung in [] (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
k = 1/(pi*e*AR);    % Faktor für Berechnung des Luftwiderstandsbeiwertes
m = 276800;         % Gewicht des Flugzeugs in [kg]
g = 9.81;           % Erdbeschleunigung in [N/kg]
alpha = 1.247015;   % Parameter zur Berechung der Luftdichte
beta = 0.000104;

% Anfangs-/End-/Zustands-/Steuerbedingungen
Z0 = [   0;         % x_0 in [m]
         0;         % h_0 in [m]
         0;         % v_0 in [m/s]
      0.27];        % gamma_0 in [rad]  
ZT = [10668;        % h_t in [m]
          0];       % gamma_t  in [rad]
q_max = 44154;      % Maximaler Staudruck in [N/m^2]
T_max = 1260000;    % Maximale Schubkraft in [N]
T_span = [0,T_max]; % Bereich des Möglichen Schub
C_L_min = 0.0;      % Minimaler Auftriebsbeiwert in []
C_L_max = 1.48;     % Maximaler Auftriebsbeiwert in []
C_L_span = [C_L_min,C_L_max]; % Bereich des möglichen Auftriebsbeiwertes
t_f = 1800;         % Endzeitpunkt in [s]
t_span = [0,t_f];   % Zeitraum

%% Problemformulierung
X = dyn_model(t_f, Z0, {T_max, C_L_min, F, C_D_0, alpha, beta, k, m, g})

Z = bvp(t_f, [Z0; Z0], {T_max, C_L_min, F, C_D_0, alpha, beta, k, m, g})


