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

% Funktionen
rho =@(h) alpha * exp(- beta * h);      % Luftdichte
q =@(v,h) 0.5 * rho(h) * v^2;           % Staudruck
C_D =@(C_L) C_D_0 + k * C_L^2;          % Luftwiderstandsbeiwert
L =@(v,h,C_L) F * C_L * q(v,h);         % Auftriebskraft
D =@(v,h,C_L) F * C_D(C_L) * q(v,h);    % Luftwiderstand

%% Problemformulierung
Z_dot =@(t,Z,T,C_L) [                                 Z(1) * cos(Z(4)); % x_dot
                                                      Z(1) * sin(Z(4)); % h_dot
                        (1/m) * (T - D(Z(3),Z(2),C_L) - m*g*sin(Z(4))); % v_dot
                     (1/(m*Z(3))) * (L(Z(3),Z(2),C_L) - m*g*cos(Z(4)))];% gamma_dot


newZ = Z_dot(1,Z0,T_max,C_L_max)
