% MaximalRangeFlight.m:
% Description:
%   Class of the optimal control problem for an Airbus A380-800, and 
%   functions for the Matlab function fmincon.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Goetz, Felix

classdef MaximalRangeFlight
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
        n_x = 5;            % Size of the state vector
        n_u = 2;            % Control vector size
        n_s = 1;            % Size of the pure state constraint vector
        n_psi = 10;          % Size of the boundary condition vector

        %% Variables for the calculation
        N;                  % Number of discretisations
        t;                  % Variable for the time vector
        ode_method;         % Implicit or explicit iterative methods to approximate the solution of a ODE
        z_0;                % Array for the start values for fmincon
        lb;                 % lower bounds for fmincon
        ub;                 % upper bounds for fmincon
    end
    methods
        %% Constructor: Set initial conditions of the direct solver
        function obj = MaximalRangeFlight(N,t,z_0,X_0,X_T,params,lb,ub,ode_method)
            % Grid size
            obj.N = N;
            obj.t = t;
            % Method for solving the ODE
            obj.ode_method = ode_method;
            % Start vector fmincon
            if obj.N == size(z_0,1)
                obj.z_0 = z_0;
            else
                obj.z_0 = z_0 .* ones(obj.N,obj.n_x+obj.n_u);
            end
            % Box conditions
            obj.lb = lb .* ones(obj.N,obj.n_x+obj.n_u);
            obj.ub = ub .* ones(obj.N,obj.n_x+obj.n_u);
            % Parameters
            obj.t_0 = params(1);
            obj.t_f = params(2);
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
        
        %% Differential equation
        function X = f(obj,t,z)
            X = [                                                                                                                 z(5)*z(4)*sind(z(2)),...
                                              z(5)*(1/(2*obj.m*z(4))) * (obj.F*z(7)*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*cosd(z(2))),...
                                                                                                                                  z(5)*z(4)*cosd(z(2)),...
                 z(5)*(1/(2*obj.m)) * (2*z(6) + (-obj.C_D_0 - obj.k*(z(7))^2)*obj.F*(z(4))^2*obj.alpha*exp(-obj.beta*z(1)) - 2*obj.m*obj.g*sind(z(2))),...
                                                                                                                                                     0];
        end
        
        %% Objective function
        function T = F_sol(obj,z)
            T = z(1,5);
        end
        
        %% Nonlinear constraint functions
        function [c,ceq] = nonlcon(obj,z)
            % Equality constraints
            c = 0.5 * obj.alpha * exp(-obj.beta*z(:,1)) .* (z(:,4)).^2 - obj.q_max;
            % Inequality constraints
            x = obj.ode_method(@obj.f,obj.t,z,obj.N,obj.n_x);
            ceq = zeros(1,obj.n_x*(obj.N-1)+obj.n_psi);
            ceq(1:obj.n_x*(obj.N-1)) = (z(1:end-1,1:obj.n_x) + x - z(2:end,1:obj.n_x))';
            ceq((obj.n_x*(obj.N-1)+1):(obj.n_x*(obj.N-1)+obj.n_psi)) = [  z(1,1)-obj.X_0(1),...
                                                                          z(1,2)-obj.X_0(2),...
                                                                          z(1,3)-obj.X_0(3),...
                                                                          z(1,4)-obj.X_0(4),...
                                                                                          0,...
                                                                        z(end,1)-obj.X_T(1),...
                                                                        z(end,2)-obj.X_T(2),...
                                                                                          0,...
                                                                                          0,...
                                                                                          0];
        end
    end
end