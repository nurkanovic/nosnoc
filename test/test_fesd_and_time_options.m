function [results,stats,model,problem_options, solver_options] = test_fesd_and_time_options(use_fesd, time_optimal_problem,...
    equidistant_control_grid, use_speed_of_time_variables, local_speed_of_time_variable)
%TEST_SIMPLE_CAR_MODEL Test the simple car model accross cross
%complementarity and mpcc modes
import casadi.*
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
solver_options.print_level = 3;
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;

% set the cross complimentarity mode that 
problem_options.cross_comp_mode = 3;
problem_options.n_s = 2;
solver_options.homotopy_update_rule = 'superlinear';
solver_options.N_homotopy = 8;

fprintf('use_fesd\ttime opt\teqdist. grid\t use sot\tlocal sot\n')
fprintf('%d\t\t\t%d\t\t\t%d\t\t\t\t\t%d\t\t\t%d\n',use_fesd, time_optimal_problem, equidistant_control_grid, use_speed_of_time_variables, local_speed_of_time_variable);
% settings being tested
problem_options.use_fesd = use_fesd;
problem_options.time_optimal_problem = time_optimal_problem;
problem_options.equidistant_control_grid = equidistant_control_grid;
problem_options.use_speed_of_time_variables = use_speed_of_time_variables;
problem_options.local_speed_of_time_variable = local_speed_of_time_variable;

% Model - define all problem functions and
% Discretization parameters
model = NosnocModel();
problem_options.N_stages = 10; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control interval (optionally a vector might be passed)
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
if time_optimal_problem
    problem_options.T = 1;    
else
    problem_options.T = 25;   
    model.f_q = u^2;
end
% Solve OCP
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
if ~isempty(results.T) && results.T < 1e-2
    warning('Something went wrong.')
    disp(results.T)
end
% disp(results.f_opt)
end

