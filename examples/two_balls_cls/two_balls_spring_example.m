clear all;
clc;
import casadi.*
close all
%%
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
% problem_options.irk_representation = 'differential';
problem_options.n_s = 3;
solver_options.print_level = 3;
% solver_options.N_homotopy = 8;
problem_options.cross_comp_mode = 1;
problem_options.dcs_mode = DcsMode.CLS;
solver_options.multiple_solvers = 0;
solver_options.mpcc_mode = "Scholtes_ineq";
% solver_options.elasticity_mode = ElasticityMode.ELL_INF;
% solver_options.psi_fun_type = CFunctionType.SCHOLTES_TWO_SIDED;
% solver_options.relaxation_method = RelaxationMode.TWO_SIDED;
problem_options.no_initial_impacts = 1;
solver_options.print_details_if_infeasible = 0;
solver_options.pause_homotopy_solver_if_infeasible = 0;
% solver_options.opts_ipopt.ipopt.linear_solver = 'ma97';
solver_options.sigma_0 = 5;
solver_options.homotopy_update_slope = 0.1;
solver_options.real_time_plot = 0;
solver_options.comp_tol  = 1e-13;
solver_options.sigma_N = 1e-13;

solver_options.mpcc_mode = "elastic_ineq";
solver_options.elastic_scholtes = 1;
solver_options.sigma_0 = 1e0;
solver_options.homotopy_update_slope = 0.2;
%%
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model = NosnocModel();
model.M = eye(2);
model.x = [q;v];
model.e = 0.8;
model.mu_f = 0;
x0 = [1;2;0;0];
model.x0 = x0;
model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
model.f_c = q(1)-R;

%% Simulation settings
N_FE = 2;
T_sim = 1;
N_sim = 200;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;

%% MATLAB solution
solver_options.use_previous_solution_as_initial_guess = 1;
[t_grid_matlab, x_traj_matlab, n_bounces, lambda_normal_guess] = two_balls_spring_matlab(T_sim,x0,model.e,1e-13);


%% Call nosnoc Integrator
initial_guess = struct();
initial_guess.x_traj = x_traj_matlab;
initial_guess.t_grid = t_grid_matlab;
initial_guess.lambda_normal_traj = lambda_normal_guess;

% [results,stats,solver] = integrator_fesd(model, settings, [], initial_guess);
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%%
plot_two_ball_traj(results);

%% compare
error = norm(x_traj_matlab(end,:)'-results.x(:,end));
fprintf('Numerical error %2.2e \n',error);


