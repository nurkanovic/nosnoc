clear all;
% clc;
import casadi.*
close all
%% init nosnoc settings and model
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%%
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
problem_options.n_s = 1;
solver_options.print_level = 3;
problem_options.cross_comp_mode = 1;
problem_options.dcs_mode = DcsMode.CLS;
problem_options.no_initial_impacts = 1;

% elasic mode with decreasing bounds for the elstaic slacks
%solver_options.mpcc_mode = "elastic_ineq";
solver_options.elastic_scholtes = 1;

solver_options.sigma_0 = 1e0;
solver_options.homotopy_update_slope = 0.1;

%% model defintion
g = 9.81;
x0 = [0.8;0];

q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 1;
model.mu_f = 0;
model.x0 = x0;
model.f_v = -g;
model.f_c = q;
%% Simulation setings
N_FE = 2;
T_sim = 3;
N_sim = 10;

problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = N_FE;

%% Analytic solution

if model.e == 0
    t_s = sqrt(2*x0(1)/g);
    tt1 = linspace(0,t_s,100);
    tt2 = linspace(t_s,T_sim,100);
    v1 = x0(2)-g*tt1;
    q1 = x0(1)+x0(2)*tt1-g*tt1.^2/2;
    if model.e == 0
        v2 = 0*(tt2-t_s);
        q2 = 0*(tt2-t_s);
    else
        v2 = -model.e*v1(end)-g*(tt2-t_s);
        q2 = q1(end)+v2(1)*(tt2-t_s)-g*(tt2-t_s).^2/2;
    end
    t_grid = [tt1,tt2];
    q_traj= [q1,q2];
    v_traj = [v1,v2];
    x_traj = [q_traj',v_traj'];
    lambda_normal = [q1*0,q2*0+g];
    Lambda_star = v2(1)-v1(end);
else
    [t_grid,x_traj,n_bounces, lambda_normal] = bouncing_ball_matlab(T_sim,model.x0,model.e,1e-12);
    Lambda_star = max(abs(diff(x_traj(:,2))));

end

%% Call nosnoc Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();

%% read and plot results
qx = results.x(1,:);
vx = results.x(2,:);

figure
subplot(311)
plot(results.t_grid,qx);
hold on
plot(t_grid,x_traj(:,1))
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
xlim([0 results.t_grid(end)])
ylim([-0.1 x0(1)+0.1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
subplot(312)
plot(results.t_with_impulse, results.x_with_impulse(2,:));
hold on
plot(t_grid,x_traj(:,2))
plot(results.t_grid,vx,'b.','MarkerSize',6);
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
ylim([-3 3])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
if problem_options.no_initial_impacts == 1
    Lm = [];
    for ii=0:length(results.Lambda_normal)-1
        if mod(ii,N_FE-1) == 0
            Lm = [Lm,nan];
        end
        Lm = [Lm, results.Lambda_normal(ii+1)];
    end
    stem(results.t_grid,[Lm,nan])
else
    stem(results.t_grid,[results.Lambda_normal,nan])
end
hold on
yline(Lambda_star,'k--')
xlim([-0.01 results.t_grid(end)])
ylim([-0.1 max([results.Lambda_normal,Lambda_star])+1])
grid on
legend({'$\Lambda$ - numerical','$\Lambda$ - analytic'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');

fprintf('position error %2.2e \n',abs(x_traj(end,1)-qx(end)))
fprintf('velocity error %2.2e \n',abs(x_traj(end,2)-vx(end)))

