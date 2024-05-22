
clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 10;
T_sim = 0.75;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'integral';
problem_options.dcs_mode = DcsMode.CLS;;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;
problem_options.no_initial_impacts = true;


model = nosnoc.model.Cls();

g = 10;
x0 = [0.8;0];

q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 1;
model.mu = 0;
model.x0 = x0;
model.f_v = -g;
model.f_c = q;

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
% %
% figure
% plot(t_grid, x_res)
% grid on
% xlabel('$t$','Interpreter','latex')
% ylabel('$x(t)$','Interpreter','latex')
% grid on



