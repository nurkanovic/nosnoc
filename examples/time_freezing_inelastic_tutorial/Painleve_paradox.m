clear all;
clear all;
clc;
import casadi.*

%% 
J = 1; % no frictioinal impulse
J = 1/32; % frictional impulse apperas
above_ground = 0.1;
%% init nosnoc
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%% settings
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 1;
problem_options.dcs_mode = 'Step';
problem_options.pss_lift_step_functions= 1;
solver_options.opts_casadi_nlp.ipopt.max_iter = 3e2;
solver_options.print_level = 2;
solver_options.N_homotopy = 10;
problem_options.cross_comp_mode = 1;
solver_options.psi_fun_type = CFunctionType.SCHOLTES;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.stagewise_clock_constraint = 0;

%%

model.e = 0;
model.mu_f = 1;
model.dims.n_dim_contact = 2;
%% the dynamics

model.a_n = 100;
qx = SX.sym('qx',1);
qy = SX.sym('qy',1);
qtheta = SX.sym('qtheta',1);
vx = SX.sym('vx',1);
vy = SX.sym('vy',1);
omega = SX.sym('omega',1);
q = [qx;qy;qtheta];
v = [vx;vy;omega];
model.x = [q;v];
model.q = q;
model.v = v;
% constraint
m = 1; l = 1;
theta0 = pi/6;
g = 9.81;
M = diag([m,m,J]);
model.M = M;
% contact points of the rod
yc = qy-l/2*cos(qtheta);
xc = qx-l/2*sin(qtheta);
model.f_v = [0;-g;0];
model.f_c = yc;
model.J_tangent = xc.jacobian(q)';
% tangent
model.x0 = [0;l/2*cos(theta0)+above_ground;theta0 ;...
           -10;0;0];
%% Simulation settings
N_finite_elements = 3;
T_sim = 0.6;
N_sim = 40;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
problem_options.use_speed_of_time_variables = 0;
problem_options.local_speed_of_time_variable = 0;
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%%
qx = results.x(1,:);
qy = results.x(2,:);
qtheta = results.x(3,:);
vx = results.x(4,:);
vy = results.x(5,:);
omega = results.x(6,:);
t = results.x(7,:);
xc_res  = [];
yc_res  = [];
for ii = 1:length(qx)
    xc_res  = [xc_res, qx(ii)-l/2*sin(qtheta(ii))];
    yc_res  = [yc_res,qy(ii)-l/2*cos(qtheta(ii))];
end
%%
h = problem_options.h_k;
figure
for ii = 1:length(qx)
    plot([qx(ii)+l/2*sin(qtheta(ii)) xc_res(ii)],[qy(ii)+l/2*cos(qtheta(ii)) yc_res(ii)],'k','LineWidth',1.5)
    hold on
    yline(0,'r')
    xlabel('$q_x$','Interpreter','latex')
    ylabel('$q_y$','Interpreter','latex')
    axis equal
    ylim([-0.1 2])
    grid on
    pause(0.01)
    clf
end
%%
figure
plot(t,vx)
hold on
plot(t,vy)
plot(t,omega)
xlabel('$t$','Interpreter','latex')
ylabel('$v$','Interpreter','latex')
grid on
