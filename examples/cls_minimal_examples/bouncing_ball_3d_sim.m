clear all;
clear all;
clc;
import casadi.*
%% init model and settings)
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%%
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 2;
solver_options.mpcc_mode = MpccMode.Scholtes_ineq;
solver_options.print_level = 3;
solver_options.N_homotopy = 10;
problem_options.use_fesd = 1;
problem_options.dcs_mode = 'CLS';
problem_options.friction_model = "Polyhedral";
problem_options.conic_model_switch_handling = "Abs";
problem_options.pss_lift_step_functions= 0;
problem_options.impose_terminal_phyisical_time  = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.nonsmooth_switching_fun = 0;
problem_options.pss_lift_step_functions = 0;
problem_options.cross_comp_mode = 1;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',3); 
v = SX.sym('v',3);
model.e = 0;
model.mu_f = 0.2;
model.dims.n_dim_contact = 3;
model.x = [q;v]; 
model.a_n = g;
model.x0 = [0;0;1;2;1;0]; 
F_ext = [1;1]*0;
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
model.D_tangent = [1,-1,0,0;
                   0,0,1,-1;
                   0,0,0,0];
%% Simulation settings
N_finite_elements = 10;
T_sim = 2;
N_sim = 1;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 0;
%% Call FESD Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%%
qx = results.x(1,:);
qy = results.x(2,:);
qz = results.x(3,:);
vx = results.x(4,:);
vy = results.x(5,:);
vz = results.x(6,:);
t_grid = results.t_grid;
figure
plot3(qx,qy,qz);
axis equal
xlim([-0.1 2])
ylim([-0.1 2])
zlim([-0.1 1])
grid on
xlabel('$q_x$','Interpreter','latex');
ylabel('$q_y$','Interpreter','latex');
zlabel('$q_z$','Interpreter','latex');
%
figure
plot(t_grid,vx);
grid on
hold on
plot(t_grid,vy);
plot(t_grid ,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$t_1^\top v$','$t_2^\top v$','$n^\top v$'},'Interpreter','latex','Location','best');
