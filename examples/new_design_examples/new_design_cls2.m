clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 10;
T_sim = 1.7;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'integral';
problem_options.dcs_mode = DcsMode.CLS;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;
problem_options.no_initial_impacts = true;
problem_options.cross_comp_mode = CrossCompMode.FE_FE;
problem_options.friction_model = "Conic";
%problem_options.friction_model = "Polyhedral";
problem_options.conic_model_switch_handling = "Abs";
%problem_options.eps_cls = 0.5;
solver_options.ipopt_callback = @cls_callback;
solver_options.sigma_0 = 10;
%solver_options.homotopy_steering_strategy = 'ELL_INF';
%solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.complementarity_tol = 1e-6;
solver_options.N_homotopy = 7;


%% model
model = nosnoc.model.Cls();

g = 10;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model.M = diag([1,1]);
model.x = [q;v];
model.e = 0;
model.mu = 0.2;
% impact
model.x0 = [0;1;4;0];
% Free flight
%model.x0 = [0;5.5;4;0];
% Friction slip slip
% model.x0 = [0;0.0;-0.6;0];
% friction slip stick
% model.x0 = [0;0.0;+0.6;0];
model.f_v = [0;-g];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.D_tangent = [1,-1;0,0];

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
% %
figure
plot(x_res(1,:), x_res(2,:))
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on



