function [results,stats,model,problem_options] = test_simple_car_model(cross_comp, mpcc_mode)
%TEST_SIMPLE_CAR_MODEL Test the simple car model accross cross
%complementarity and mpcc modes
import casadi.*
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
% set the cross complimentarity mode that 
problem_options.cross_comp_mode = cross_comp;
solver_options.mpcc_mode = mpcc_mode;
fprintf('cross_comp\tmpcc_mode\n')
fprintf('%d\t\t%s\n', cross_comp, mpcc_mode);
problem_options.n_s = 2;
solver_options.N_homotopy = 10; 
% Time-settings  - Solve an time optimal control problem
problem_options.time_optimal_problem = 1;
% Model - define all problem functions and
% Discretization parameters
model = NosnocModel();
problem_options.N_stages = 10; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;    % Time horizon
% Symbolic variables and bounds
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value
% bounds on states
model.lbx = [-inf;-20];
model.ubx = [inf;20];
% control
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];
% Solve OCP
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
end

