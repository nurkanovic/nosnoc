clear all
close all
import casadi.*
import vdx.*

R = 3.5;
%% Define projected system
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
%problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
%problem_options.rk_scheme = RKSchemes.LOBATTO_IIIC;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
%problem_options.rk_representation= RKRepresentation.differential_lift_x; 
problem_options.rk_representation = RKRepresentation.integral;
problem_options.cross_comp_mode = CrossCompMode.STAGE_STAGE;
%problem_options.cross_comp_mode = CrossCompMode.FE_FE;
problem_options.N_finite_elements = 2;
problem_options.n_s = 2;
problem_options.N_stages = 25;
problem_options.T = 1;
problem_options.rho_h = 1e-4;
problem_options.time_optimal_problem = true;
problem_options.use_fesd = true;
problem_options.gcs_lift_gap_functions = false;
%solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.complementarity_tol = 1e-10;
solver_options.print_level = 3;

model = nosnoc.model.Pds();
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
x = [x1;x2];
x_target = [0;0;0;0];
model.x = [x];
model.lbx = [-inf;-inf;-inf;-inf];
model.ubx = [inf;inf;inf;inf];
x0 =[-25;-25;-15;-15];
model.x0 = [x0];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
model.u = [u1;u2];
model.lbu = [-100/sqrt(2);-100/sqrt(2);-60/sqrt(2);-60/sqrt(2)];
model.ubu = [100/sqrt(2);100/sqrt(2);60/sqrt(2);60/sqrt(2)];
%model.lbu = [-100/sqrt(2);-100/sqrt(2);0;0];
%model.ubu = [100/sqrt(2);100/sqrt(2);0;0];
model.u0 = model.ubu;
model.c = [norm_2(x2-x1)-2*R];
model.f_x_unconstrained = [u1;u2];

% costs
model.f_q = 0;
model.f_q_T = 0.5*(x-x_target)'*(x-x_target);

% Solve
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% plot
fontsize = 12;
x_res = ocp_solver.get('x');
u_res = ocp_solver.get('u');
u_rep = kron(u_res, ones(1,problem_options.N_finite_elements(1)));
T_final = ocp_solver.get('T_final');
lambda_res = ocp_solver.get('lambda');
t_res = ocp_solver.get_time_grid();
if ~problem_options.use_fesd
    t_res = t_res*T_final
end
c_fun = casadi.Function('nabla_c', {x}, {model.c});
nabla_c_fun = casadi.Function('nabla_c', {x}, {model.c.jacobian(x)'});
c_res = full(c_fun(x_res(1:4,:)));
nabla_c_res = full(nabla_c_fun(x_res(1:4,:)));
%plot_discs(h_res,x_res,[3.5,3.5], ["circle", "circle"])
v_res = u_rep + repmat(lambda_res(2:end),4,1).*nabla_c_res(:, 2:end);
figure('Position', [0,0, 400. 400])
ax = subplot(3,1,1);
plot(t_res(1:end-1), v_res([1,3],:), "LineWidth", 2)
xlabel("$t$")
ylabel("$v_x$")
ylim([40, 75])
xlim([0, T_final])
grid on
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ax = subplot(3,1,2);
plot(t_res(1:end-1), v_res([2,4],:), "LineWidth", 2)
xlabel("$t$")
ylabel("$v_y$")
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ylim([40, 75])
xlim([0, T_final])
grid on
ax = subplot(3,1,3);
plot(t_res, c_res, "LineWidth", 2)
xlabel("$t$")
ylabel("$c(x)$")
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ylim([-0.1, 10])
xlim([0, T_final])
grid on
