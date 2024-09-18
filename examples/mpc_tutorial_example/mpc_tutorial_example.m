clear; clc; close all;
import casadi.*
% 
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.solver.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation. 

% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)

% Time-settings  - Solve a 
problem_options.N_stages = 6; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 3; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 3;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)

% define differential states and populate the model.
q = SX.sym('q'); 
v = SX.sym('v'); % CasADi symbolic variables for states
model.x = [q;v]; % populate model state vectors
model.x0 = [0;0]; % initial value
v_max = 20;
model.lbx = [-inf;-v_max]; % lower bounds on states
model.ubx = [inf;v_max]; % upper bounds on states
% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
model.u = u;
u_max = 5;
model.lbu = -u_max ; 
model.ubu = u_max ;
% Dynamics of the piecewise smooth systems
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
% Define the regions of the PSS
v_threshold = 10;
model.c = v-v_threshold; % single switching functions (implies two regions)
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [f_1 f_2]; % The columns of this matrix store the vector fields of every region.

model.f_q = (q-200)^2; % Add terminal quadratic cost
model.f_q_T = (q-200)^2 + 0.1*(v)^2; % Add terminal quadratic cost

% Setup solver an mpc options
solver_options.homotopy_update_rule = 'superlinear'; % Use superlinear update rule for relaxation parameter sigma.
solver_options.homotopy_update_slope = 0.05; % Rate which the relaxation sigma is reduced at: sigma_i+1 = kappa*sigma_i 
solver_options.homotopy_update_exponent = 2; % Rate which the relaxation sigma is reduced at: sigma_i+1 = sigma_i^kappa
solver_options.complementarity_tol = 1e-7; % Value to drive the complementarity residual to.
solver_options.N_homotopy = 10; % Maximum number of homotopy iterations.

% Set the fast sigma for followup solves to the complementarity tol to only do 1 nlp solve. 
mpc_options.fullmpcc_fast_sigma_0 = 1e-7;

% create mpc object
mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, solver_options);

% Do MPC assuming the predicted state is accurate, in practice this may not be true.
x = model.x0; u = []; t = 0;
x0 = x;
for step=1:50
    [u_i, stats] = mpc.get_feedback(x0);
    mpc.do_preparation();
    x0 = mpc.get_predicted_state();
    x = [x, x0];
    u = [u,u_i];
    t = [t, t(end) + problem_options.h];
end

% Plot
figure
subplot(311)
plot(t,x(1,:))
xlabel("$t$")
ylabel("$q$")
ylim([-5 205])
subplot(312)
plot(t,x(2,:))
xlabel("$t$")
ylabel("$v$")
ylim([-5 25])
subplot(313)
stairs(t,[u,u(end)])
xlabel("$t$")
ylabel("$u$")
ylim([-5.5 5.5])
