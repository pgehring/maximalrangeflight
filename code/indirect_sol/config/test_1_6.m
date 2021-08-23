% test_1_6.m

% Versuchsaufbau:
%   - Veränderte Endzeit

%% Speicher Parameter
results_name = 'test_1_6';

%% Testparameter + Lösungsmethode der ODE und Objekt der Problemklasse erhalten
N = 5;

direkt_sol = readmatrix(strcat('test_1_1','.txt'))';
eta_mehrsv = zeros(8,N);
for k = 0:N-1
    eta_mehrsv(:,k+1) = [direkt_sol(1:4,k+1);5.10402542449782;1145.48666826056;-1;-2.60344266831659];
end

eta_mehrsv = [-0.000461697366532953,0.986114562416963,3.94764693225948,7.71668739410971,11.0625721620951;-0.0846487572569101,0.271832230846204,0.869521184399733,0.842654452877248,0.976943535547555;4.40587449745183e-17,307.290839190894,620.813036849306,919.679653888999,1212.15142003352;99.9978160790220,104.680947420473,101.064854423120,98.5519710076505,96.4818863792967;-2.65200000524707,-2.74740160969817,-2.79123668783286,-2.79226931489557,-2.79323154274863;-1383.47152361742,-891.400583555400,-15.8162930540916,959.031785341558,1924.40581959055;-0.999999999999994,-0.999999999999996,-0.999999999999997,-0.999999999999998,-1.00000000000000;-18.9352714810131,-4.61674473792164,3.74850152514068,5.65595239889774,4.45341656466738];

eta_mehrsv = [-0.000481975506538138,0.983633607838976,3.94109842143356,7.71249326544199,11.0609796116599;-0.100385357970859,0.243353535948076,0.835620487551300,0.811961408861671,0.957998534032198;-3.14355978083808e-17,307.290326478921,620.816087567980,919.692173411492,1212.17274665505;99.9976264622472,104.680970517861,101.067283327747,98.5551759262831,96.4847792014170;-2.65180370311138,-2.74719096635784,-2.79101259997609,-2.79204457122570,-2.79300673345979;-1383.29644371719,-891.225409245816,-15.6836652166738,959.099329408668,1924.40581959055;-1.00000000000001,-1.00000000000000,-1.00000000000000,-0.999999999999999,-0.999999999999998;-18.9323109955578,-4.61548360945551,3.74798860098606,5.65517027738106,4.45285364513099];

eta_mehrsv = [-0.000485285018069210,0.983229435948301,3.94003509821920,7.71181980676347,11.0607232420631;-0.102954781310113,0.238703550850367,0.830085096512906,0.806949722427508,0.954905108662823;-2.83348569960466e-17,307.290244892354,620.816601469418,919.694253482993,1212.17628376616;99.9975955208845,104.680974300159,101.067683634756,98.5557022874376,96.4852555339296;-2.65177151270314,-2.74715643101093,-2.79097587182775,-2.79200773626022,-2.79296988778308;-1383.26779651703,-891.196814161113,-15.6620490171351,959.110323251314,1924.40581959055;-0.999999999999986,-0.999999999999985,-0.999999999999984,-0.999999999999987,-0.999999999999992;-18.9318282842141,-4.61527879090751,3.74790405932726,5.65504211447884,4.45276130207080];


X_0 = [   0;           % h_0 in [m]
       0.27;           % gamma_0 in [rad]  (Steigflug mit Neigungswinkel von cs 20°)
          0;           % x_0 in [m]
        100];          % v_0 in [m/s] (Benötigte Startgeschwindigkeit zum Abheben)

X_T = [10668;          % h_t in [m]
           0];         % gamma_t  in [Grad]
X_T = direkt_sol(1:2,N);
       
params = [       0,... % t_0:   Anfangszeitpunkt in [s]
               N*3,... % t_f:   Endzeitpunkt in [s]
          1.247015,... % alpha: Parameter zur Berechung der Luftdichte in []
          0.000104,... % beta: 
              9.81,... % g:     Erdbeschleunigung in [N/s^2]
             0.032,... % C_D_0: Nullluftwiderstandsbeiwert in []             (Der Luftwiderstandsbeiwert ohne Auftrieb ist ein dimensionsloser Parameter, der die Luftwiderstandskraft eines Flugzeugs mit seiner Größe, Geschwindigkeit und Flughöhe in Beziehung setzt)
               0.8,... % e:     Oswaldfaktor in []                           (Der Oswald-Wirkungsgrad ist ein Korrekturfaktor, der die Änderung des Luftwiderstands mit dem Auftrieb eines dreidimensionalen Flügels oder Flugzeugs im Vergleich zu einem idealen Flügel mit demselben Seitenverhältnis darstellt)
               845,... % F:     Wirksame Fläche in [m^2]                     (von der Luft angeströmte Fläche) 
               7.5,... % AR:    Flügelstreckung in []                        (aspect ratio) (Das Seitenverhältnis eines Flügels ist definiert als das Verhältnis seiner Spannweite zu seiner mittleren Sehne) -> (Spannweite in [m])^2 / Tragflächeninhalt in [m^2] = b^2 / F
            276800,... % m:     Leergewicht des A380 in [kg]
             44154];   % q_max: Maximaler Staudruck in [N/m^2]  

prob = MaximalRangeFlightIndirect(eta_mehrsv,X_0,X_T,params);

%% Boxbeschränkungen
AbsTol        =  1e-8;
RelTol        =  1e-8;
StopTol       = 1e-12;
StopTolArmijo = 1e-15;
maxit = 100;

% ode_method = @ode45;
ode_method = @ode23s;
shooting_methods = shooting_methods(ode_method,0,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,' ');
shooting_method = @shooting_methods.Mehrfachschiessverfahren;