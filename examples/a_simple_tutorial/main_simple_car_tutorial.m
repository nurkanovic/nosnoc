clear all
import casadi.*
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.cross_comp_mode = 7;
problem_options.n_s = 2;
% Time-settings  - Solve an time optimal control problem
problem_options.time_optimal_problem = 1;
problem_options.N_stages = 10; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;    % Time horizon

% settings.nlpsol = 'snopt';  % Note: requires installing.

% Model - define all problem functions and
% Discretization parameters
model = nosnoc.model.Pss();
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
% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% Plot
v_max = 20;
u_max = 5;
v_trash_hold = 10;
plot_results_nosnoc_tutorial
