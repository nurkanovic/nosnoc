close all
clear all
import casadi.*

problem_options = NosnocProblemOptions();
problem_options.n_s = 2;
problem_options.dcs_mode = 'Step';
problem_options.print_level = 7;

model = NosnocModel();
x1 = SX.sym('x1');
model.x = x1;
model.c = [x1];
model.S = [-1; 1];

f_11 = 3;
f_12 = 1;
model.F = [f_11, f_12];
model.x0 = -1;

model.T = pi/4;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 2;
problem_options.cross_comp_mode = 3;
model.dims.n_s = 2;

mpcc = NosnocMPCC(problem_options, model.dims, model);

solver_options = NosnocSolverOptions();
solver_options.mpcc_mode = MpccMode.ell_1_penalty;

solver = NosnocSolver(mpcc, solver_options);

%solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
