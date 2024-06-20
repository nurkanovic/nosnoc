clear all;
clc;
import casadi.*
close all
%% init
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
solver_options.print_level = 3;
solver_options.N_homotopy = 7;
problem_options.cross_comp_mode = 7;
problem_options.dcs_mode = DcsMode.CLS;
problem_options.no_initial_impacts = 1;
solver_options.complementarity_tol = 1e-6;
%problem_options.friction_model = "Polyhedral";
problem_options.friction_model = "Conic"; % "Conic"
problem_options.conic_model_switch_handling = "Lp";
% settings.mpcc_mode = MpccMode.elastic_ineq;
% solver_options.nlpsol = 'snopt';
solver_options.use_previous_solution_as_initial_guess = 1;
%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model.M = diag([1,1]);
model.x = [q;v];
model.e = 0;
model.mu_f = 0.2;
% impact
model.x0 = [0;1;4;0];
% Free flight
% model.x0 = [0;5.5;4;0];
% Friction slip slip
% model.x0 = [0;0.0;-0.6;0];
% friction slip stick
% model.x0 = [0;0.0;+0.6;0];
model.f_v = [3*1;-g];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.D_tangent = [1,-1;0,0];

%% Simulation settings
N_FE = 2;
T_sim = 1.7;
N_sim = 10;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;

%% Call nosnoc Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%% read and plot results
unfold_struct(results,'base');
qx = results.x(1,:);
qy = results.x(2,:);
vx = results.x(3,:);
vy = results.x(4,:);
figure(1)
subplot(121)
plot(qx,qy);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(results.t_grid,vy);
hold on
plot(results.t_grid,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
%% forces
if N_sim == 1
figure
subplot(311)
plot(results.t_grid,results.x(3:4,:));
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(312)
plot(results.t_grid,[nan*ones(model.n_contacts,1),results.lambda_normal])
hold on
plot(results.t_grid,[nan*ones(model.n_tangents,1),results.lambda_tangent])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda$','interpreter','latex');
subplot(313)
stem(results.t_grid,[results.Lambda_normal,nan*ones(model.n_contacts,1)])
hold on
stem(results.t_grid,[results.Lambda_tangent,nan*ones(model.n_tangents,1)]')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');
xlim([0 T_sim])
end
