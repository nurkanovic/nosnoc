clear all;
clear all;
clc;
import casadi.*
close all
%%
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
problem_options.dcs_mode = DcsMode.Step;
problem_options.n_s = 3;
problem_options.cross_comp_mode = 1;
solver_options.print_level = 2;
problem_options.use_fesd = 1;
solver_options.comp_tol = 1e-5;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e4;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 0;

problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 0;
problem_options.stagewise_clock_constraint = 0;

%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',1); 
v = SX.sym('v',1);
model = NosnocModel();
model.x = [q;v]; 
model.e = 0.9;
model.k_aux = 50;
model.x0 = [1;0]; 
model.f_v = -g;
model.f_c = q;
model.dims.n_dim_contact = 1;
%% Simulation settings
N_FE = 3;
T_sim = 3;
N_sim = 60;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%% read and plot results
unfold_struct(results,'base');
qx = results.x(1,:);
vx = results.x(2,:);
t_opt = results.x(3,:);
figure
subplot(121)
plot(t_opt,qx);
axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$q$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vx);
hold on
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

%% speed of time
figure
subplot(121)
plot(results.t_grid,t_opt)
hold on
plot(results.t_grid,results.t_grid,'k--')
grid on
xlabel('$\tau$','interpreter','latex');
ylabel('$t$','interpreter','latex');
subplot(122)
stairs(results.s_sot)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');

%% complementarity residuals
