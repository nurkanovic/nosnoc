%clear all;
clc;
import casadi.*
close all
%% init nosnoc
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%%
N_FE = 3;
T_sim = 0.2;
N_sim = 1;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.cross_comp_mode = 1;
%problem_options.lift_complementarities = 1;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.pss_lift_step_functions = 0;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.n_s = 2;


%%
solver_options.complementarity_tol = 1e-10;
solver_options.sigma_N = 1e-10;
solver_options.N_homotopy = 16;
solver_options.print_level = 3;
solver_options.store_integrator_step_results = 1;
solver_options.use_previous_solution_as_initial_guess = 0;
%solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%%
g = 10;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); 
v = SX.sym('v',2);

model.x = [q;v]; 
model.e = 0.0;
model.mu = 0.0;
model.x0 = [0;0.06;3;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.dims.n_dim_contact = 2;

model.verify_and_backfill(problem_options);

tf_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, problem_options)
