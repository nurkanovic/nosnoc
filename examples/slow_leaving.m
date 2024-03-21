close all
clear all
import casadi.*

problem_options = NosnocProblemOptions();
problem_options.dcs_mode = 'Stewart';
problem_options.print_level = 1;

model = NosnocModel();
x = SX.sym('x', 2);
model.x = x;
model.c = [x(2)+0.25];
model.S = [-1; 1];

f = [x(2);-x(1)];
f_11 = f + model.c.jacobian(x)';
f_12 = f;
model.F = [f_11, f_12];
model.x0 = [1.0;-.25];

problem_options.T_sim = 5.0;
problem_options.N_sim = 27;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 2;
problem_options.cross_comp_mode = 7;
problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;
problem_options.n_s = 2;

solver_options = NosnocSolverOptions();
solver_options.comp_tol = 1e-12;
solver_options.N_homotopy = 15;
solver_options.sigma_N = 1e-15;
solver_options.mpcc_mode = MpccMode.elastic_ineq;
solver_options.opts_casadi_nlp.ipopt.max_iter = 20000;

%mpcc = NosnocMPCC(problem_options, model);
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();

plot(results.t_grid,results.x)
grid on
