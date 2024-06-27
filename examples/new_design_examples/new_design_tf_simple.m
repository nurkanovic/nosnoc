%clear all;
clc;
import casadi.*
close all
%% init nosnoc
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%%
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.dcs_mode = DcsMode.Heaviside;
problem_options.cross_comp_mode = 7;
%problem_options.lift_complementarities = 1;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.s_sot_max = 25;
%problem_options.rho_sot = 1e-1;
problem_options.rho_h = 1;
problem_options.k_aux = 100;
problem_options.relax_terminal_physical_time = ConstraintRelaxationMode.NONE;
problem_options.rho_terminal_physical_time = 1e4;
%problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.NONE;
%problem_options.rho_terminal_numerical_time = 1e3;
%problem_options.step_equilibration = StepEquilibrationMode.linear_complementarity;

solver_options.print_level = 3;
solver_options.complementarity_tol = 1e-10;
solver_options.sigma_N = 1e-10;
solver_options.N_homotopy = 16;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_opts.use_previous_solution_as_initial_guess = 0;

%%
g = 10;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); 
v = SX.sym('v',2);

model.x = [q;v]; 
model.e = 0;
model.mu = 0.0;
model.x0 = [0;1.0;1;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.dims.n_dim_contact = 2;

%% Simulation settings
N_FE = 4;
T_sim = 1.1;
N_sim = 21;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.n_s = 3;

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
% %
figure
plot(x_res(end,:)', x_res(1:2,:)')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
figure
plot(x_res_full(end,:)', x_res_full(1:2,:)')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
figure
plot(x_res(end,:), x_res(3:4,:))
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
grid on
figure
plot(x_res_full(1,:), x_res_full(2,:))
figure
plot(x_res_full(2,:), x_res_full(4,:))
figure
hold on
plot(t_grid, x_res(end,:))
plot(t_grid, t_grid, 'k--')
hold off
