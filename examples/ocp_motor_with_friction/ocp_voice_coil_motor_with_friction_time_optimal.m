%% Info
% This is an optimal control example from the paper:
% Optimal control of a voice-coil-motor with Coulombic friction
% by Bahne Christiansen; Helmut Maurer; Oliver Zirn
% Published in: 2008 47th IEEE Conference on Decision and Control
% DOI: 10.1109/CDC.2008.4739025
clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % 
solver_options = nosnoc.solver.Options();

%% set some options
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.N_stages = 30; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 0.09;    % Time horizon
problem_options.cross_comp_mode = 'FE_FE';
problem.options.dcs_mode = "Heaviside"; % or "Heaviside"
problem_options.time_optimal_problem = true;
problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;
problem_options.use_fesd = 1;

% problem_options.gamma_h = 0;
% problem_options.N_finite_elements = 1;
% problem_options.use_fesd = 1;
% problem_options.print_level  = 5;
% problem_options.equidistant_control_grid = 0;

% solver options
solver_options.N_homotopy = 7;
solver_options.homotopy_update_rule = 'superlinear';
solver_options.opts_casadi_nlp.ipopt.max_iter = 1000;
% solver_options.print_level =  5;
%% Create model
% model = nosnoc.model.stewart();
model = nosnoc.model.Pss(); 
% Parameters
m1 = 1.03; % slide mass
m2 = 0.56; % load mass
k = 2.4e3; %  spring constant N/m
c = 0.00; % damping
U_max = 5; % voltage Back-EMF, U = K_s*v_1;
R  = 2; % coil resistancel ohm
L = 2e-3; % inductivity, henry
K_F = 12; % force constant N/A ; F_L = K_F*I; % Lorenz force
K_S = 12; % Vs/m (not provided in the paper above)
F_R = 2.1; % guide friction forec, N
x0 = [0;0;0;0;0];
x_target = [0.01;0;0.01;0;0];
x1 = SX.sym('x1');  % motor mass position
v1 = SX.sym('v1'); % motor mass velocity
x2 = SX.sym('x2');  % load mass position
v2 = SX.sym('v2'); % load mass velocity
I = SX.sym('I'); % electrin current
x = [x1;v1;x2;v2;I];
u = SX.sym('U'); % the motor voltage
model.x = x;
model.x0 =  x0;         

model.u = u;
model.lbu = -U_max;
model.ubu = U_max;

A = [0  1   0   0   0;...
    -k/m1 -c/m1 k/m1 c/m1 K_F/m1;...
    0   0   0   1   0;...
    k/m2   c/m2   -k/m2   -c/m2 0;...
    0 -K_S/L    0       0   -R/L];
B = [zeros(4,1);1/L];
C1 = [0;-F_R/m1;0;0;0]; % v1 >0
C2 = -C1; %v1<0
% switching dynamics with different friction froces
f_1 = A*x+B*u+C1; % v1>0
f_2 = A*x+B*u+C2; % v1<0
F = [f_1, f_2];
c1 = v1;  % Switching function
model.S = [1;-1];
model.F = F;
model.c = c1;
model.g_terminal = x-x_target;

%% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%% Plot results
x = ocp_solver.get('x');
u = ocp_solver.get('u');
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();

x1_opt = x(1,:);
v1_opt= x(2,:);
x2_opt= x(3,:);
v2_opt= x(4,:);
I_opt= x(5,:);

figure
subplot(411)
plot(t_grid,x1_opt)
hold on
plot(t_grid,x2_opt)
ylabel('$x(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$x_1(t)$','$x_2(t)$'},'Interpreter','latex','Location','best')
subplot(412)
plot(t_grid,v1_opt)
hold on
plot(t_grid,v2_opt)
yline(0,'k:');
ylabel('$v(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$v_1(t)$','$v_2(t)$'},'Interpreter','latex','Location','best')
subplot(413)
plot(t_grid,I_opt)
ylabel('$I(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
% t_grid_u = t_grid_u';
subplot(414)
stairs(t_grid_u,[u,nan]);
ylabel('$u(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
ylim([-6.1 6.1])
