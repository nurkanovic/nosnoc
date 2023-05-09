% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

%
%
% TODO: remove varargout / varargin and make API clear!
% Examples of calling the function
% [results] = integrator_fesd(model,settings);
% [results,stats] = integrator_fesd(model,settings);
% [results,stats,model] = integrator_fesd(model,settings);

function [varargout] = integrator_fesd(model, settings, u_sim)

control_exists = 1;
if ~exist('u_sim','var')
    control_exists = 0;
    u_sim = [];
end


results.t_grid = [];
results.x_res = [];
stats = [];

%% generate time-freezing model before turning off time-related settings
if settings.time_freezing && ~solver_exists
    [model,settings] = time_freezing_reformulation(model,settings);
end

%% Settings of integration
[model] = refine_model_integrator(model,settings);

%% Create solver functions for integrator step
solver = NosnocSolver(model, settings);
model = solver.model;
settings = solver.settings;

% TODO remove this but needed now because unfold_struct clobers fields for some reason
% TODO when moving to model class remove this
settings_bkp = settings;
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
settings = settings_bkp;

%% check does the provided u_sim has correct dims
if control_exists
    [n_rows,n_cols] = size(u_sim);
    if n_rows~=n_u || n_cols~=N_sim
        error('Matrix with control inputs has the wrong size. Required dimension is n_u x N_sim.')
    end
end
%% Initialization
x_res = [x0];
x_res_extended = [x0];
diff_res = [];
alg_res = [];
% Stewart
theta_res = [];
lambda_res = [];
mu_res = [];
% Output including all RK-stage values
lambda_res_extended = [];
theta_res_extended = [];
mu_res_extended = [];

% Step
alpha_res = [];
lambda_0_res = [];
lambda_1_res = [];
% Output including all RK-stage values
alpha_res_extended = [];
lambda_0_res_extended = [];
lambda_1_res_extended = [];

% step size
h_vec = [];
% sot (in time-freezing)
s_sot_res = [];
% states
complementarity_stats  = [];
homotopy_iteration_stats = [];
time_per_iter = [];
simulation_time_pased = 0;
W = [];
all_res = [];
%% Main simulation loop
for ii = 1:N_sim
    if control_exists
        % NOTE: This should work because sim only has N_stages=0
        solver.set('u', u_sim(:,ii));
    end

    [sol,stats] = solver.solve();
    res = extract_results_from_solver(model, solver.problem, settings, sol);
    all_res = [all_res,res];
    time_per_iter = [time_per_iter; stats.cpu_time_total];
    % verbose
    simulation_time_pased = simulation_time_pased + model.T;
    if print_level >=2
        fprintf('Integration step %d / %d (%2.3f s / %2.3f s) converged in %2.3f s. \n',ii, N_sim,simulation_time_pased,T_sim,time_per_iter(end));
    end
    % Store differential states
    w_opt = full(sol.nlp_results(end).x);
    W = [W, w_opt];
    diff_states = w_opt(solver.problem.ind_x_all);
    alg_states = w_opt(solver.problem.ind_z_all);

    diff_res = [diff_res;diff_states];
    alg_res = [alg_res;alg_states];

    % step-size
    h_opt = w_opt(flatten_ind(solver.problem.ind_h));
    % differential
    x_opt_extended = w_opt(solver.problem.ind_x_all);
    x_opt_extended  = reshape(x_opt_extended,n_x,length(x_opt_extended)/n_x);

    % only bounadry value
    if isequal(irk_representation,'integral')
        x_opt  = res.x_opt(:,2:end);
    elseif isequal(irk_representation, 'differential_lift_x')
        x_opt  = res.x_opt(:,2:end);
    else
        x_opt  = res.x_opt(:,2:end);
    end


    % TODO: this should use indices instead of n_*
    switch dcs_mode
        case 'Stewart'
            alg_states_extended = reshape(alg_states,n_z_all,length(alg_states)/n_z_all);
            theta_opt_extended = [alg_states_extended(1:n_theta,:)];
            lambda_opt_extended = [alg_states_extended(n_theta+1:2*n_theta,:)];
            mu_opt_extended = [alg_states_extended(end-n_sys+1:end,:)];

            theta_opt= theta_opt_extended(:,1:n_s:end);
            lambda_opt= lambda_opt_extended(:,1:n_s:end);
            mu_opt= mu_opt_extended(:,1:n_s:end);
        case 'Step'
            alg_states_extended = reshape(alg_states,n_z_all,length(alg_states)/n_z_all);
            alpha_opt_extended = [alg_states_extended(1:n_alpha,:)];
            lambda_0_opt_extended = [alg_states_extended(n_alpha+1:2*n_alpha,:)];
            lambda_1_opt_extended = [alg_states_extended(2*n_alpha+1:3*n_alpha,:)];

            alpha_opt= alpha_opt_extended(:,1:n_s:end);
            lambda_0_opt= lambda_0_opt_extended(:,1:n_s:end);
            lambda_1_opt= lambda_1_opt_extended(:,1:n_s:end);
    end

    % update initial guess and inital value
    x0 = x_opt(:,end);
    %     update clock state
    if impose_terminal_phyisical_time
        solver.problem.p0(end) = solver.problem.p0(end)+model.T;
    end
    solver.set("x0", x0);
    

    % TODO Set up homotopy solver to take p_val explicitly
    if use_previous_solution_as_initial_guess
        % TODO make this possible via solver interface directly
        solver.problem.w0(n_x+1:end) = w_opt(n_x+1:end);
    end

    % Store data
    if use_fesd
        h_vec = [h_vec;h_opt];
    else
        h_vec = [h_vec;h_k(1)*ones(N_stages*N_finite_elements(1),1)];
    end
    %sot
    s_sot_res  = [s_sot_res,w_opt(flatten_ind(solver.problem.ind_sot))];
    %differntial.
    x_res = [x_res, x_opt(:,end-N_finite_elements(1)*N_stages+1:end)];
    x_res_extended = [x_res_extended,x_opt_extended(:,2:end)];

    % algebraic
    switch dcs_mode
      case 'Stewart'
        theta_res = [theta_res, theta_opt];
        lambda_res = [lambda_res, lambda_opt];
        mu_res = [mu_res, mu_opt];
        theta_res_extended = [theta_res_extended,theta_opt_extended ];
        lambda_res_extended = [lambda_res_extended,lambda_opt_extended];
        mu_res_extended = [mu_res_extended,mu_opt_extended];

      case 'Step'
        alpha_res = [alpha_res, alpha_opt];
        lambda_0_res = [lambda_0_res, lambda_0_opt];
        lambda_1_res = [lambda_1_res, lambda_1_opt];

        alpha_res_extended = [alpha_res_extended, alpha_opt_extended];
        lambda_0_res_extended = [lambda_0_res_extended, lambda_0_opt_extended];
        lambda_1_res_extended = [lambda_1_res_extended, lambda_1_opt_extended];
    end
    %stats
    complementarity_stats  = [complementarity_stats; stats.complementarity_stats(end)];
    homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];

    % plot during execution
    if real_time_plot
        if time_freezing
            figure(100)
            clf
            plot(x_res(end,:),x_res(1:end-1,:));
            xlabel('$t$ [phyisical time]','Interpreter','latex');
            ylabel('$x(t)$','Interpreter','latex');
            grid on
            xlim([0 T_sim]);
            ylim([min(min(x_res(1:end-1,:)))-0.3 max(max(x_res(1:end-1,:)))+0.3])
            hold on
        else
            figure(100)
            clf
            t_temp = [0,cumsum(h_vec)'];
            plot(t_temp,x_res(:,1:end));
            xlabel('$t$','Interpreter','latex');
            ylabel('$x(t)$','Interpreter','latex');
            grid on
            xlim([0 T_sim]);
            ylim([min(x_res(:))-0.3 max(x_res(:))+0.3])
            hold on
        end
    end
end
total_time = sum(time_per_iter);
%% Verbose
fprintf('\n');
fprintf('----------------------------------------------------------------------------------------------------------------------\n');
if use_fesd
    fprintf( ['Simulation with the FESD ' char(irk_scheme) ' with %d-RK stages completed.\n'],n_s);
else
    fprintf( ['Simulation with the standard ' char(irk_scheme) ' with %d-RK stages completed.\n'],n_s);
end
fprintf( ['RK representation: ' char(irk_representation) '.\n']);
fprintf('---------------------------------------------- Stats summary----------------------------------------------------------\n');
fprintf('N_sim\t step-size\t\tN_stg\tN_FE\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter\tMax. comp.\tMin. comp.\n');
fprintf('%d\t\t\t%2.3f\t\t%d\t\t%d\t\t%2.3f\t\t\t\t%2.3f\t\t\t%2.3f\t\t\t\t%2.2e\t%2.2e\n',N_sim,h_sim,N_stages,N_finite_elements(1),total_time,max(time_per_iter),min(time_per_iter),max(complementarity_stats),min(complementarity_stats));
fprintf('----------------------------------------------------------------------------------------------------------------------\n\n');
%% Output

results.x_res  = x_res;
% Output all stage values as well.
results.x_res_extended  = x_res_extended;
switch dcs_mode
  case 'Stewart'
    results.theta_res = theta_res;
    results.lambda_res = lambda_res;
    results.mu_res = mu_res;
    results.theta_res_extended  = theta_res_extended;
    results.lambda_res_extended  = lambda_res_extended;
    results.mu_res_extended  = mu_res_extended;
  case 'Step'
    results.alpha_res = alpha_res;
    results.lambda_0_res = lambda_0_res;
    results.lambda_1_res = lambda_1_res;
    results.alpha_res_extended = alpha_res_extended;
    results.lambda_0_res_extended = lambda_0_res_extended;
    results.lambda_1_res_extended = lambda_1_res_extended;
  case 'CLS'
    % TODO
end
stats.complementarity_stats = complementarity_stats;
stats.time_per_iter = time_per_iter;
stats.homotopy_iteration_stats = homotopy_iteration_stats;

results.h_vec = h_vec;
results.s_sot_res = s_sot_res;
results.t_grid = cumsum([0;h_vec])';

% full diff and alg variables (e.g. for initalizing a solver)
results.diff_res = diff_res;
results.alg_res = alg_res;
results.W = W;
results.solver_ouput = sol;
results.all_res = all_res;

varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
varargout{5} = solver;
end

