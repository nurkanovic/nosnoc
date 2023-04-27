clear all;
clear all;
clc;
import casadi.*

%% 
J = 1; % no frictioinal impulse
J = 1/32; % frictional impulse apperas
above_ground = 0.1;
%%
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 3;
settings.dcs_mode = 'CLS';
settings.friction_model = "Polyhedral";
settings.conic_model_switch_handling = "Abs";
settings.pss_lift_step_functions= 0;
settings.opts_ipopt.ipopt.max_iter = 3e2;
settings.print_level = 3;
settings.N_homotopy = 5;
settings.cross_comp_mode = 1;
settings.time_freezing = 0;
settings.impose_terminal_phyisical_time = 1;
settings.mpcc_mode = "elastic_ineq";
%%
model.e = 0;
model.mu = 0;
model.n_dim_contact = 2;
%% the dynamics
model.n_q = 3;
model.g = -9.81*1;
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
yc_left = qy-l/2*cos(qtheta);
xc_left = qx-l/2*sin(qtheta);
yc_right = qy+l/2*cos(qtheta);
xc_right = qx+l/2*sin(qtheta);
model.f_v = [0;-g;0];
model.f_c = [yc_left;yc_right];
model.J_tangent = xc_left.jacobian(q)';
model.D_tangent = [xc_left.jacobian(q)',-xc_left.jacobian(q)'];

% tangent
model.x0 = [0;l/2*cos(theta0)+above_ground;theta0 ;...
           -10;0;0];

above_ground = 0.18;
theta0 = 0.75*pi/2;
model.x0 = [0;l/2*cos(theta0)+above_ground;theta0 ;...
           0;0;0];

%% Simulation setings
N_finite_elements = 50;
T_sim = 0.8;
N_sim = 1;
model.T_sim = T_sim;
model.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call FESD Integrator
if settings.time_freezing
    [model,settings] = time_freezing_reformulation(model,settings);
    settings.time_freezing = 0;
end
[results,stats,model] = integrator_fesd(model,settings);
%%
qx = results.x_res(1,:);
qy = results.x_res(2,:);
qtheta = results.x_res(3,:);
vx = results.x_res(4,:);
vy = results.x_res(5,:);
omega = results.x_res(6,:);
if settings.time_freezing
    t = results.x_res(7,:);
else
    t = results.t_grid;
end

xc_res  = [];
yc_res  = [];
for ii = 1:length(qx)
    xc_res  = [xc_res, qx(ii)-l/2*sin(qtheta(ii))];
    yc_res  = [yc_res,qy(ii)-l/2*cos(qtheta(ii))];
end
%%
h = model.h_k;
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
    pause(h)
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
