clear all;
clear all;
clc;
import casadi.*
%% init nosnoc
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();  
model = NosnocModel();
%% settings
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 1;
solver_options.mpcc_mode = MpccMode.Scholtes_ineq;
solver_options.print_level = 3;
solver_options.N_homotopy = 6;
problem_options.use_fesd = 1;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions= 1;
problem_options.impose_terminal_phyisical_time  = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.nonsmooth_switching_fun = 1;
problem_options.pss_lift_step_functions = 0;
solver_options.use_previous_solution_as_initial_guess = 0;
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
% norm(F_ext)
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
%% Simulation settings
N_finite_elements = 3;
T_sim = 3;
N_sim = 20;
problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = N_finite_elements;

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
t_opt = results.x(7,:);
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
%%
figure
plot(t_opt,vx,'LineWidth',2);
grid on
hold on
plot(t_opt,vy);
plot(t_opt,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$t_1^\top v$','$t_2^\top v$','$n^\top v$'},'Interpreter','latex','Location','best');


%%
t_grid = results.t_grid;
lambda0 = results.lambda_n;
lambda1 = results.lambda_p;
alpha = results.alpha;

if problem_options.time_freezing_nonlinear_friction_cone
theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:));
theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:));
theta = [theta1;theta2;theta3];
else
theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(1-alpha(4,:));
theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(alpha(4,:));
theta4  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(1-alpha(4,:));
theta5  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(alpha(4,:));
theta = [theta1;theta2;theta3;theta4;theta5];
end

dims = model.dims
n_f = dims.n_theta_step;
t_grid(1) = [];
figure
for ii = 1:n_f 
    subplot(1,n_f,ii)
    plot(t_grid,theta(ii,:));
    grid on
    xlabel('$\tau$','Interpreter','latex')
    ylabel(['$\theta_' num2str(ii) '$'],'Interpreter','latex')
    ylim([-0.1 1.1])
end
%% switching functions
% c_fun = model.c_fun;
% x_res = results.x;
% t_grid = results.t_grid;
% c_eval = [];
% for ii = 1:length(x_res)
%     c_eval = [c_eval,full(c_fun(x_res(:,ii)))];
% end
% figure
% plot(t_grid,c_eval))
% grid on

