close all
clear all
import casadi.*

problem_options = NosnocProblemOptions();
problem_options.n_s = 1;
problem_options.dcs_mode = 'Stewart';
problem_options.print_level = 1;

model = NosnocModel();
x1 = SX.sym('x1');
model.x = x1;
model.c = [x1];
model.S = [-1; 1];

f_11 = 3;
f_12 = 1;
model.F = [f_11, f_12];
model.x0 = -1;

problem_options.T = pi/4;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 2;
problem_options.cross_comp_mode = 1;
%problem_options.lift_complementarities = 1;
problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;
problem_options.n_s = 2;

solver_options = nosnoc.solver.Options();
%solver_options.mpcc_mode = MpccMode.elastic_ineq;
%solver_options.psi_fun_type = CFunctionType.STEFFENSON_ULBRICH;
%solver_options.elasticity_mode = ElasticityMode.ELL_INF;
solver_options.opts_casadi_nlp.ipopt.max_iter = 20000;

mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();

plot(results.t_grid,results.x)
grid on
