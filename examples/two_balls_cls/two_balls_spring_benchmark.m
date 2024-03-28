clear all;
clc;
import casadi.*
close all

%%
benchmark_globals;

%% create model
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);


%% create reference MATLAB solution
ref_sol_filename = "two_balls_guess_sol.mat";
[t_grid_guess, x_traj_guess, n_bounces_guess, lambda_normal_guess] = two_balls_spring_matlab(1.1*T_sim, x0, e, 1e-3);
save(ref_sol_filename, "t_grid_guess", "x_traj_guess", "n_bounces_guess", "lambda_normal_guess");
% load(ref_sol_filename)
tic

%% run experiments
for irk_scheme = IRK_SCHEMES
    for with_guess = [0]
        for n_s = NS_VALUES
            for N_sim = NSIM_VALUES
                for N_FE = NFE_VALUES
                    model = NosnocModel();
                    model.M = eye(2);
                    model.x = [q;v];
                    model.e = e;
                    model.mu_f = 0;
                    model.x0 = x0;
                    model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
                    model.f_c = q(1)-R;
                    % settings
                    problem_options = NosnocProblemOptions();
                    solver_options = nosnoc.solver.Options();
                    problem_options.irk_scheme = irk_scheme;
                    % problem_options.irk_representation = 'differential';
                    problem_options.n_s = n_s;
                    solver_options.print_level = 3;
                    problem_options.cross_comp_mode = 1;
                    problem_options.dcs_mode = DcsMode.CLS;
                    solver_options.multiple_solvers = 0;
                    problem_options.no_initial_impacts = 1;
                    solver_options.print_details_if_infeasible = 0;
                    solver_options.pause_homotopy_solver_if_infeasible = 0;
                    % solver_options.opts_ipopt.ipopt.linear_solver = 'ma97';
                    solver_options.comp_tol  = 1e-13;
                    solver_options.sigma_N = 1e-13;
                    solver_options.mpcc_mode = "elastic_ineq";
                    solver_options.elastic_scholtes = 1;
                    solver_options.sigma_0 = 1e0;
                    solver_options.homotopy_update_slope = 0.2;
                    solver_options.opts_casadi_nlp.ipopt.max_iter = 1500;
                    solver_options.use_previous_solution_as_initial_guess  = 1;
                    %% Simulation settings
                    problem_options.T_sim = T_sim;
                    problem_options.N_finite_elements = N_FE;
                    problem_options.N_sim = N_sim;

                    %% Call nosnoc Integrator
                    initial_guess = struct();
                    initial_guess.x_traj = x_traj_guess;
                    initial_guess.t_grid = t_grid_guess;
                    initial_guess.lambda_normal_traj = lambda_normal_guess;

                    if with_guess
                        % settings.opts_casadi_nlp.ipopt.least_square_init_duals = 'yes';
                        integrator = NosnocIntegrator(model, problem_options, solver_options, [], initial_guess);
                        [results,stats] = integrator.solve();
                    else
                        integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
                        [results,stats] = integrator.solve();
                    end

                    results_filename = get_results_filename(n_s, N_sim, N_FE, problem_options.irk_scheme, with_guess);
                    save(results_filename, 'results', 'stats', 'problem_options', 'solver_options')

                    clear model solver
                end
            end
        end
    end
end
disp('experiment loop took')
toc


%% read and plot results
q1 = results.x(1,:);
q2 = results.x(2,:);
v1 = results.x(3,:);
v2 = results.x(4,:);
t_grid = results.t_grid;

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
ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
Lambda_opt = [results.Lambda_normal];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan])
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');

%% compare
error = norm(x_traj_guess(end,:)'-results.x(:,end));
fprintf('Numerical error %2.2e \n',error);
