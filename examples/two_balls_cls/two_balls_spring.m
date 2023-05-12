clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
% settings.irk_representation = 'differential';
settings.n_s = 3;
settings.print_level = 3;
% settings.N_homotopy = 8;
settings.cross_comp_mode = 3;
settings.dcs_mode = DcsMode.CLS;
settings.multiple_solvers = 0;
settings.mpcc_mode = "Scholtes_ineq";
settings.no_initial_impacts = 1;
% settings.opts_ipopt.ipopt.linear_solver = 'ma97';
settings.sigma_0 = 5;
settings.homotopy_update_slope = 0.1;

%%
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model.M = eye(2);
model.x = [q;v];
model.e = 0.8;
model.mu = 0;
x0 = [1;2;0;0];
model.x0 = x0;
model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
model.f_c = q(1)-R;

%% Simulation settings
N_FE = 2;
T_sim = 1;
N_sim = 1000;

model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;

%% MATLAB solution
[t_grid_matlab, x_traj_matlab, n_bounces] = two_springs_matlab(T_sim,x0,model.e,1e-5);


%% Call nosnoc Integrator
initial_guess = struct();
initial_guess.x_traj = x_traj_matlab;
initial_guess.t_grid = t_grid_matlab;
settings.sigma_0 = 1e-3;

[results,stats,model,settings,solver] = integrator_fesd(model, settings, [], initial_guess);

%% read and plot results
unfold_struct(results,'base');
q1 = x_res(1,:);
q2 = x_res(2,:);
v1 = x_res(3,:);
v2 = x_res(4,:);

%%
figure
subplot(311)
plot(t_grid,q1,'LineWidth',1.5);
hold on
plot(t_grid,q2,'LineWidth',1.5);
yline(R,'k--')
xlim([0 t_grid(end)])
% ylim([-1.0 max([q1,q2])+1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(312)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
xlim([0 t_grid(end)])
% ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
Lambda_opt = [results.all_res.Lambda_normal_opt];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan])
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');

%% compare
error = norm(x_traj_matlab(end,:)'-x_res(:,end));
fprintf('Numerical error %2.2e \n',error);


