% MaximalRangeFlight.m:
% Description:
%   Class of the optimal control problem for an Airbus A380-800, and 
%   functions for the Matlab function fmincon.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Goetz, Felix

classdef maximal_range_flight
    properties
        %% Parameters for the optimal control problem
        t_0                 % Starting time in [s]
        t_f                 % End time in [s]

        X_0 = [   0;        % h_0 in [m]
               0.27;        % gamma_0 in [rad]
                  0;        % x_0 in [m]
                100];       % v_0 in [m/s]

        X_T = [10668;       % h_t in [m]
                   0];      % gamma_t in [Grad]

        alpha = 1.247015;   % Parameters for calculating the air density in []
        beta = 0.000104;    % Parameters for calculating the air density in []
        g = 9.81;           % Gravitational acceleration in [N/s^2]
        C_D_0 = 0.032;      % Zero air resistance coefficient in []             
        e = 0.8;            % Oswaldfactor in []
        F = 845;            % Effective area in [m^2]
        AR = 7.5;           % Aspect ratio in []
        k;                  % Factor for calculation of the drag coefficient in []
        m = 276800;         % Empty weight of the A380 in [kg]
        q_max = 44154;      % Maximum dynamic pressure in [N/m^2]
        T_min = 0;          % Minimum thrust in [N]
        T_max = 1260000;    % Maximum thrust in [N]
        C_L_min = 0.0;      % Minimum lift coefficient in []
        C_L_max = 1.48;     % Maximum lift coefficient in []

        %% Problem sizes
        n_x = 4;            % Size of the state vector
        n_u = 2;            % Control vector size
        n_s = 1;            % Size of the pure state constraint vector
        n_psi = 8;          % Size of the boundary condition vector

        %% Variables for the calculation
        N;                  % Number of discretisations
        tspan;                  % Variable for the time vector
        z_0;                % Array for the start values for fmincon
    end
    methods
        %% Constructor: Set initial conditions of the direct solver
        function obj = maximal_range_flight(z_0,X_0,X_T,params)
            % Start vector fmincon
            obj.z_0 = z_0;
            % Parameters
            obj.t_0 = params(1);
            obj.t_f = params(2);
            obj.tspan = [obj.t_0,obj.t_f];
            obj.X_0 = X_0;
            obj.X_T = X_T;
            obj.alpha = params(3);
            obj.beta = params(4);
            obj.g = params(5);
            obj.C_D_0 = params(6);
            obj.e = params(7);
            obj.F = params(8);
            obj.AR = params(9);
            obj.k = 1/(pi*obj.e*obj.AR);
            obj.m = params(10);
            obj.q_max = params(11);
        end

        %% Funktionen Randwertproblem
        function [T,C_L] = control(obj,t,z)
            % Toleranz für = 0 Bedingung 
            tol = 1e-16;
            % Synthese-Steuerung: Auftriebsbeiwert
            K_1 = (obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(6))/(2*obj.m);
            K_2 = (obj.k*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
            %
            b_1_1 = (K_1<0 && K_2<0 && K_1/(2*K_2)<obj.C_L_max);
            b_1_2 = (K_1<0 && K_2<0 && K_1/(2*K_2)>obj.C_L_max);
            b_2 = (abs(K_1)<tol && K_2<0);
            b_3 = (K_1>0 && K_2<0);
            %
            b_4 = (K_1<0 && abs(K_2)<tol);
            % b_5 = (abs(K_1)<tol && abs(K_2)<tol);
            b_6 = (K_1>0 && abs(K_2)<tol);
            %
            b_7 = (K_1<0 && K_2>0);
            b_8 = (abs(K_1)<tol && K_2>0);
            b_9_1 = (K_1>0 && K_2>0 && (K_1/K_2)<obj.C_L_max);
            b_9_2 = (K_1>0 && K_2>0 && (K_1/K_2)>obj.C_L_max);
            b_9_3 = (K_1>0 && K_2>0 && abs((K_1/K_2)-obj.C_L_max)<tol);
            %
            if b_2 || b_3 || b_6 || b_9_2
                C_L = obj.C_L_min;
            elseif b_1_1
                C_L = K_1/(2*K_2);
            elseif b_1_2 || b_4 || b_7 || b_8 || b_9_1
                C_L = obj.C_L_max;
            elseif b_9_3
                C_L = obj.C_L_max;
            else
                C_L = (obj.C_L_max-obj.C_L_min)/2;
            end
            % Synthese-Steuerung: Schub
            if z(8) > 0
                T = obj.T_min;
            elseif z(8) < 0
                T = obj.T_max;
            else
                T = (obj.T_max-obj.T_min)/2;
            end
        end
        
        function Z_dot = G(obj,t,z)   
%             % Lösung von C_L mit fmincon 
% %             H =@(u) (obj.F*obj.alpha*exp(-obj.beta*z(1))*u(2)*z(4)*z(6))/(2*obj.m) + (u(1)*z(8))/obj.m - ((obj.C_D_0+obj.k*u(2)^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
% %             H =@(u) (obj.F*obj.alpha*exp(-obj.beta*z(1))*u(2)*z(4)*z(6))/(2*obj.m) + (u(1)*z(8))/obj.m - (obj.k*u(2)^2*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
% %             H =@(u) K_1*u(2) + (u(1)*z(8))/obj.m - K_2*u(2)^2;
%             H =@(u) sind(z(2))*z(4)*z(5) +...
%                     (obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(6))/(2*obj.m)*u(2) +...
%                     -(obj.g*cosd(z(2))*z(6))/(z(4)) ...
%                     + cosd(z(2))*z(4)*z(7)...
%                     +(u(1)*z(8))/obj.m ...
%                     - (obj.k*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m)*u(2)^2 ...
%                     - (obj.C_D_0*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m)...
%                     - obj.g*sind(z(2))*z(8);
% 
%             start = [0;1.48];
% %             start = [T;C_L];
%             lb = [obj.T_min;obj.C_L_min];
%             ub = [obj.T_max,obj.C_L_max];
% %             options = optimoptions('fmincon','Display','off','ConstraintTolerance',1e-8,'StepTolerance',1e-11);
%             options = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxFunctionEvaluations',2000.0e+03,'MaxIterations',4.0e+05,'ConstraintTolerance',1e-8,'StepTolerance',1e-11,'OptimalityTolerance',1e-10);
%             [u_sol,fval,exitflag,output] = fmincon(H,start,[],[],[],[],lb,ub,[],options);
%             
% %             T = u_sol(1);
% %             C_L = u_sol(2);
%             
%             if abs(T-u_sol(1)) > 1
%                 diff = T-u_sol(1)
%                 H([T;C_L])
%                 H(u_sol)
%                 disp('ERROR')
%             elseif abs(C_L-u_sol(2)) > 1e-3
%                 diff = C_L - u_sol(2)
%                 disp('ERROR 2')
%             elseif (exitflag ~= 1) && (abs(H([T;C_L])-H(u_sol))>10)
%                 H([T;C_L])
%                 H(u_sol)
%                 disp('ERROR 3')
%             end
%             
%             plot = 0;
%             if plot ==1
%                 diff_T = linspace(obj.T_min,obj.T_max,50);
%                 diff_C_L = linspace(obj.C_L_min,obj.C_L_max,50);
%                 n_T = length(diff_T);
%                 n_C_L = length(diff_C_L);
%                 sol_grid = zeros(n_T,n_C_L);
%                 for i=1:length(diff_T)
%                     for j=1:length(diff_C_L)
%                         sol_grid(j,i) = H([diff_T(i),diff_C_L(j)]);
%                     end
%                 end                
%                 surf(diff_T,diff_C_L,sol_grid)
%             end
            
            % Get control
            [T,C_L] = control(obj,t,z);
            % 
            h_dot = z(4)*sind(z(2));
            gamma_dot = (obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L*z(4))/(2*obj.m) - (obj.g*cosd(z(2)))/z(4);
            x_dot = z(4)*cosd(z(2));
            v_dot = T/obj.m - ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2)/(2*obj.m) - obj.g*sind(z(2));
            
            lambda1_dot = + (obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(4)*z(6))/(2*obj.m)...
                            - ((obj.C_D_0+obj.k*C_L^2)*obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
            lambda2_dot = - cosd(z(2))*z(4)*z(5) - (obj.g*sind(z(2))*z(6))/z(4)...
                            + sind(z(2))*z(4)*z(7) + cosd(z(2))*obj.g*z(8);
            lambda3_dot = 0;
            lambda4_dot = - sind(z(2))*z(5) - ((obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L)/(2*obj.m)+(obj.g*cosd(z(2)))/(z(4)^2))*z(6)...
                            - cosd(z(2))*z(7) + ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(8))/obj.m;
            %
            Z_dot = [h_dot;gamma_dot;x_dot;v_dot;lambda1_dot;lambda2_dot;lambda3_dot;lambda4_dot;];    
        end

        function Z = G_Z(obj,t,z)
            % Synthese-Steuerung
            [~,C_L] = control(obj,t,z);
            %
            J_1_2 = z(4)*cosd(z(2));
            J_1_4 = sind(z(2));
            %
            J_2_1 = -(obj.F*obj.alpha*obj.beta*exp(-obj.beta*z(1))*z(4)*C_L)/(2*obj.m);
            J_2_2 = (obj.g*sind(z(2)))/z(4);
            J_2_4 = (obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L)/(2*obj.m) + (obj.g*cosd(z(2)))/(z(4)^2);
            %
            J_3_2 = -z(4)*sind(z(2));
            J_3_4 = cosd(z(2));
            %
            J_4_1 = ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*obj.beta*exp(-obj.beta*z(1))*z(4)^2)/(2*obj.m);
            J_4_2 = -obj.g*cosd(z(2));
            J_4_4 = -((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4))/(obj.m);
            %
            J_5_1 = (obj.alpha*obj.beta^2*obj.F*exp(-obj.beta*z(1))*C_L*z(4)*z(6))/(2*obj.m) - ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.beta^2*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
            J_5_4 = -(obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(6))/(2*obj.m) + ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.beta*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(8))/(obj.m);
            J_5_6 = -(obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(4))/(2*obj.m);
            J_5_8 = ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.beta*obj.alpha*exp(-obj.beta*z(1))*z(4)^2)/(2*obj.m);
            %
            J_6_2 = -sind(z(2))*z(4)*z(5) + (obj.g*cosd(z(2))*z(6))/z(4) - cosd(z(2))*z(4)*z(7) + sind(z(2))*obj.g*z(8);
            J_6_4 = cosd(z(2))*z(5) - (obj.g*sind(z(2))*z(6))/(z(4)^2) - sind(z(2))*z(7);
            J_6_5 = cos(z(2))*z(4);
            J_6_6 = (obj.g*sind(z(2)))/z(4);
            J_6_7 = -sind(z(2))*z(4);
            J_6_8 = -cosd(z(2))*obj.g;
            %
            J_8_1 = -(obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(6))/(2*obj.m) + ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.beta*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(8))/(obj.m);
            J_8_2 = cosd(z(2))*z(5) - (obj.g*sind(z(2))*z(6))/(z(4)^2) - sind(z(2))*z(7);
            J_8_4 = -(2*obj.g*cosd(z(2))*z(6))/(z(4)^3) - ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(8))/(obj.m);
            J_8_5 = sind(z(2));
            J_8_6 = (C_L*obj.F*obj.alpha*exp(-obj.beta*z(1)))/(2*obj.m) + (obj.g*cosd(z(2)))/(z(4)^2);
            J_8_7 = cosd(z(2));
            J_8_8 = -((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4))/(obj.m);
            %
            Z = [    0,J_1_2,    0,J_1_4,    0,    0,    0,    0;
                 J_2_1,J_2_2,    0,J_2_4,    0,    0,    0,    0;
                     0,J_3_2,    0,J_3_4,    0,    0,    0,    0;
                 J_4_1,J_4_2,    0,J_4_4,    0,    0,    0,    0;
                 J_5_1,    0,    0,J_5_4,    0,J_5_6,    0,J_5_8;
                     0,J_6_2,    0,J_6_4,J_6_5,J_6_6,J_6_7,J_6_8;
                     0,    0,    0,    0,    0,    0,    0,    0;
                 J_8_1,J_8_2,    0,J_8_4,J_8_5,J_8_6,J_8_7,J_8_8];
        end

        function z = R(obj,z_0,z_f)
            lambda_0 = 1; % Auszugehen von lambda_0 = 0 oder lambda_0 = 1
            % 
            z = [z_0(1)-obj.X_0(1);
                 z_0(2)-obj.X_0(2);
                 z_0(3)-obj.X_0(3);
                 z_0(4)-obj.X_0(4);
                 z_f(1)-obj.X_T(1);
                 z_f(2)-obj.X_T(2);
                   z_f(7)+lambda_0;
                          z_f(8)-0];
            
%             Konstanz der Hamilton-Funktion
%             [T_0,C_L_0] = control(obj,1,z_0);
%             [T_f,C_L_f] = control(obj,1,z_f);
%             H =@(u,z) sind(z(2))*z(4)*z(5) +...
%                     (obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(6))/(2*obj.m)*u(2) +...
%                     -(obj.g*cosd(z(2))*z(6))/(z(4)) ...
%                     + cosd(z(2))*z(4)*z(7)...
%                     +(u(1)*z(8))/obj.m ...
%                     - (obj.k*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m)*u(2)^2 ...
%                     - (obj.C_D_0*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m)...
%                     - obj.g*sind(z(2))*z(8);
%             diff = H([T_0,C_L_0],z_0) - H([T_f,C_L_f],z_f);
% 
%             z = [z_0(1)-obj.X_0(1);
%                  z_0(2)-obj.X_0(2);
%                  z_0(3)-obj.X_0(3);
%                  z_0(4)-obj.X_0(4);
%                  z_f(1)-obj.X_T(1);
%                  z_f(2)-obj.X_T(2);
%                    z_f(7)+lambda_0;
%                           z_f(8)-0;
%                               diff];
        end

        function z = R_Z_0(obj,z_0,z_f)
            z = [1,0,0,0,0,0,0,0;
                 0,1,0,0,0,0,0,0;
                 0,0,1,0,0,0,0,0;
                 0,0,0,1,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0];
            
%             z=z_0;
%             [T,C_L] = control(obj,1,z_0);
%             
%             H_h_dot = + (obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(4)*z(6))/(2*obj.m)...
%                             - ((obj.C_D_0+obj.k*C_L^2)*obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
%             H_gamma_dot = - cosd(z(2))*z(4)*z(5) - (obj.g*sind(z(2))*z(6))/z(4)...
%                             + sind(z(2))*z(4)*z(7) + cosd(z(2))*obj.g*z(8);
%             H_x_dot = 0;
%             H_v_dot = - sind(z(2))*z(5) - ((obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L)/(2*obj.m)+(obj.g*cosd(z(2)))/(z(4)^2))*z(6)...
%                             - cosd(z(2))*z(7) + ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(8))/obj.m;
%             H_lambda1_dot = sind(z(2))*z(4);
%             H_lambda2_dot = (obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L*z(4))/(2*obj.m) - (obj.g*sind(z(2)))/z(4);
%             H_lambda3_dot = cosd(z(2))*z(4);
%             H_lambda4_dot = T/obj.m - ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2)/(2*obj.m) - obj.g*sind(z(2));
%             
%             z = [1,0,0,0,0,0,0,0;
%                  0,1,0,0,0,0,0,0;
%                  0,0,1,0,0,0,0,0;
%                  0,0,0,1,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  H_h_dot,H_gamma_dot,H_x_dot,H_v_dot,H_lambda1_dot,H_lambda2_dot,H_lambda3_dot,H_lambda4_dot];
        end

        function z = R_Z_f(obj,z_0,z_f)
            z = [0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0;
                 1,0,0,0,0,0,0,0;
                 0,1,0,0,0,0,0,0;
                 0,0,0,0,0,0,1,0;
                 0,0,0,0,0,0,0,1];
            
%             z=z_f;
%             [T,C_L] = control(obj,1,z_f);
%             
%             H_h_dot = + (obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*C_L*z(4)*z(6))/(2*obj.m)...
%                             - ((obj.C_D_0+obj.k*C_L^2)*obj.alpha*obj.beta*obj.F*exp(-obj.beta*z(1))*z(4)^2*z(8))/(2*obj.m);
%             H_gamma_dot = - cosd(z(2))*z(4)*z(5) - (obj.g*sind(z(2))*z(6))/z(4)...
%                             + sind(z(2))*z(4)*z(7) + cosd(z(2))*obj.g*z(8);
%             H_x_dot = 0;
%             H_v_dot = - sind(z(2))*z(5) - ((obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L)/(2*obj.m)+(obj.g*cosd(z(2)))/(z(4)^2))*z(6)...
%                             - cosd(z(2))*z(7) + ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)*z(8))/obj.m;
%             H_lambda1_dot = sind(z(2))*z(4);
%             H_lambda2_dot = (obj.F*obj.alpha*exp(-obj.beta*z(1))*C_L*z(4))/(2*obj.m) - (obj.g*sind(z(2)))/z(4);
%             H_lambda3_dot = cosd(z(2))*z(4);
%             H_lambda4_dot = T/obj.m - ((obj.C_D_0+obj.k*C_L^2)*obj.F*obj.alpha*exp(-obj.beta*z(1))*z(4)^2)/(2*obj.m) - obj.g*sind(z(2));
%             
%             z = [0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  0,0,0,0,0,0,0,0;
%                  1,0,0,0,0,0,0,0;
%                  0,1,0,0,0,0,0,0;
%                  0,0,0,0,0,0,1,0;
%                  0,0,0,0,0,0,0,1;
%                  H_h_dot,H_gamma_dot,H_x_dot,H_v_dot,H_lambda1_dot,H_lambda2_dot,H_lambda3_dot,H_lambda4_dot];
        end
        
        %% Calculating the control for the plots
        function Z = sol_func(obj,sol)
            for i = 1:length(sol)
                [T,C_L] = control(obj,1,sol(i,:));
                Z(i,:) = [sol(i,1:4),T,C_L];
            end            
        end
        
    end
end