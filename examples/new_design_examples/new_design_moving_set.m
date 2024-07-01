
clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 100;
T_sim = 10;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 3;
%problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
%problem_options.rk_representation= 'differential_lift_x'; 
problem_options.rk_representation= 'integral';
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;

%solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.print_level = 2;

model = nosnoc.model.Pds();

model.x0 = [0;pi-2;0];
x = SX.sym('x',2);
t = SX.sym('t');
model.x = [x;t];
model.c = [-norm(x - [sin(t);cos(t)])^2+(1-0.5)^2];
model.f_x_unconstrained = [0; 0; 1];
model.E = diag([1,1,0]);

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate('natural_residual_ineq');
%
figure
plot(x_res(1,:), x_res(2,:))
grid on
xlabel('$x_1(t)$','Interpreter','latex')
ylabel('$x_2(t)$','Interpreter','latex')
grid on

fig = figure('Position', [10 10 1600 800]);
plot_moving_set([],x_res,[0.5], ["circle"], fig, 'moving_set');
