clear all
close all
import casadi.*
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
% Choosing the Runge - Kutta Method and number of stages
problem_options.n_s = 2; % number of runge Kutta stages
problem_options.dcs_mode = "Stewart";
model = nosnoc.model.Pss();
N_stages = 50; % number of data points 
N_FE = 2; % number of intermediate integration steps between data points
problem_options.N_stages = N_stages; % number of control intervals, or data poitns in this context
problem_options.N_finite_elements = N_FE; % number of integration steps (in one control intevral)
problem_options.T = 2;    % Time/simulation horizon
problem_options.step_equilibration = StepEquilibrationMode.l2_relaxed;
% Symbolic variables
x = SX.sym('x');
model.x = x;
model.x0 = -1; % inital value
% Dyanmics and the regions
model.c = x; % swiching function
f_1 = 1;  % for c < 0 , hence -1 in first entry of S
f_2 = -1;  % for c > 0 , hence 1 in second entry of S
model.S = [-1;1];
model.F = [f_1 f_2]; % collect all dynamics modes in one matrix

ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

x_res = ocp_solver.get("x");
t_grid = ocp_solver.get_time_grid();

%%
figure
plot(t_grid, x_res)
grid on
xlabel('t')
ylabel('x(t)')
hold on
%% Set up estimation problem
% create data by perturbing simulation results with random points from [-0.1,0.1]
dx = 0.1*1;
x_data = x_res + (-dx + 2*dx.*rand(size(x_res)));
plot(t_grid,x_data,'x')
x_samples = x_data(1:N_FE:end); 

%% Create nosnoc estimation problem
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
% Choosing the Runge - Kutta Method and number of stages
problem_options.n_s = 2; % number of runge Kutta stages
problem_options.dcs_mode = "Stewart";
problem_options.step_equilibration = 'direct';
problem_options.euler_cost_integration = 1;
model = nosnoc.model.Pss();
problem_options.N_stages = N_stages; % number of control intervals (no controls, set to one)
problem_options.N_finite_elements = N_FE; % number of integration steps (in one control intevral)
problem_options.T = 2;    % Time/simulation horizon
% Symbolic variables
x = SX.sym('x');
x_data_sym = SX.sym('x_data_sym'); % symbolic variable for data
v_sys = SX.sym('v_sys',2); % Uknown system parameters (optimization variables)
model.x = x;
model.v_global = v_sys; % name global, because they are not time dependent
model.ubv_global = [10;10];
model.lbv_global = [-10;-10];
model.x0 = -1; % inital value
% Dyanmics and the regions
model.c = x; % swiching function
% determine two parameters (dynamics parameteres, for the two modes of the system)
f_1 = v_sys(1);  % for c < 0 , hence -1 in first entry of S
f_2 = v_sys(2);  % for c > 0 , hence 1 in second entry of S
model.S = [-1;1];
model.F = [f_1 f_2]; % collect all dynamics modes in one matrix
model.p_time_var = x_data_sym;
model.p_time_var_val = x_samples(2:end); % value of the parameters
model.f_q = 1e3*(x-x_data_sym)'*(x-x_data_sym); % linear least squares objective

estimator = nosnoc.ocp.Solver(model, problem_options, solver_options);
estimator.solve();

x_est = estimator.get("x");
t_grid_est = estimator.get_time_grid();

%%
plot(t_grid_est, x_est)

