function [results] = polish_homotopy_solution(model,problem,settings,results,eps_sigma)
% This functions fixes heuristicly the active set from the provided initial
% guess (e.g., from the solution of the last NLP in the homotopy loop) and
% solves a smooth NLP (without complementarity conditions).
% If you call this function outside of the nosnoc_solver, please provide
% the settings struct from the nosnoc_solver model output (as it may
% contain useful information, e.g. the time-freezing model has been already
% autogenerated).
%% input
eps_sigma = 10*eps_sigma;


%% alert
if isequal(settings.dcs_mode,'Stewart');
    error('Polishing mode is currently avilable for settings.dcs_mode = ''Step''.');
end
%% unfold
import casadi.*
unfold_struct(settings,'caller');
unfold_struct(model,'caller');


%% proces raw results to prepare for active sets
% if exist('results.x_opt','var')
    % check has it been processed already
    results = extract_results_from_solver(model,problem,settings,results);
% end
unfold_struct(results,'caller');
%% thake x on relevant points, and get intex sets
% x mathaching on points where the algebraic variables are defined
x_modified = x_opt_extended;
x_modified(:,1:n_s+1:end) = [];
% getting the index sets of the algebraic variables
ind_z_extended = reshape(ind_z,n_z,length(ind_z)/n_z);
ind_alpha= [ind_z_extended(1:n_alpha,:)];
ind_lambda0 = [ind_z_extended(n_alpha+1:2*n_alpha,:)];
ind_lambda1  = [ind_z_extended(2*n_alpha+1:3*n_alpha,:)];
%% evaluating switching futction
c_eval = [];
for ii = 1:length(x_modified)
    % TODO Make this work with parameters
    c_eval = [c_eval,full(c_fun(x_modified(:,ii),[]))];
end

if polishing_derivative_test
    z_opt_extended =  w_opt(ind_z_extended);
    dot_c_eval = [];
    n_inner = N_finite_elements(1)*n_s;
    for ii = 1:N_stages
        for jj = 1:n_inner
            ii_temp = (ii-1)*n_inner+jj;
            if n_u > 0
                dot_c_eval = [dot_c_eval,full(dot_c_fun(x_modified(:,ii_temp),z_opt_extended(:,ii_temp),u_opt(:,ii)))];
            else
                dot_c_eval = [dot_c_eval,full(dot_c_fun(x_modified(:,ii_temp),z_opt_extended(:,ii_temp)))];
            end
        end
    end
end
%% Check all signs

ind_negative = c_eval < -eps_sigma;
ind_positive = c_eval > eps_sigma ;
ind_zero = abs(c_eval) <= eps_sigma;
%% Check time derivative of |c|<eps
eps_sigma_derivative = eps_sigma*10;
if polishing_derivative_test
    ind_sliding_zero= abs(dot_c_eval)<= eps_sigma_derivative;
    ind_sliding_positive  = dot_c_eval>eps_sigma_derivative;
    ind_sliding_negative = dot_c_eval<-eps_sigma_derivative;

    ind_sliding_zero = ind_zero & ind_sliding_zero;
    ind_sliding_positive  = ind_zero & ind_sliding_positive;
    ind_sliding_negative = ind_zero & ind_sliding_negative;
end

%% Asign values
if polishing_derivative_test
    % c > 0 , lambda0 = 0, alpha = 1;
    ind_lambda0_fixed = ind_lambda0(ind_positive | ind_sliding_positive);
    ind_alpha1_fixed = ind_alpha(ind_positive | ind_sliding_positive);
    % c < 0 , lambda1 = 0, alpha = 0;
    ind_lambda1_fixed = ind_lambda1(ind_negative | ind_sliding_negative);
    ind_alpha0_fixed = ind_alpha(ind_negative | ind_sliding_negative);
    % c = 0 , lambda1 = 0, lambda0 = 0;
    % her we must make clear is it sliding or is it corssing
    ind_lambda1_sliding = ind_lambda1(ind_sliding_zero);
    ind_lambda0_sliding = ind_lambda0(ind_sliding_zero);
else
    % c > 0 , lambda0 = 0, alpha = 1;
    ind_lambda0_fixed = ind_lambda0(ind_positive);
    ind_alpha1_fixed = ind_alpha(ind_positive);
    % c < 0 , lambda1 = 0, alpha = 0;
    ind_lambda1_fixed = ind_lambda1(ind_negative);
    ind_alpha0_fixed = ind_alpha(ind_negative);
    % c = 0 , lambda1 = 0, lambda0 = 0;
    % her we must make clear is it sliding or is it corssing
    ind_lambda1_sliding = ind_lambda1(ind_zero);
    ind_lambda0_sliding = ind_lambda0(ind_zero);
end

%% create nlp
settings.mpcc_mode = 'direct';
settings.opts_ipopt.ipopt.tol = 1e-14;
settings.opts_ipopt.ipopt.max_iter = 1.5e3;
solver = NosnocSolver(model, settings);

%% set appropiate bounds to zero
solver.problem.lbw(ind_lambda1_fixed(:)) = 0; solver.problem.ubw(ind_lambda1_fixed(:)) = 0;
solver.problem.lbw(ind_alpha0_fixed(:)) = 0; solver.problem.ubw(ind_alpha0_fixed(:)) = 0;

solver.problem.lbw(ind_lambda0_fixed(:)) = 0;  solver.problem.ubw(ind_lambda0_fixed(:)) = 0;
solver.problem.lbw(ind_alpha1_fixed(:)) = 1;  solver.problem.ubw(ind_alpha1_fixed(:)) = 1;
%
solver.problem.lbw(ind_lambda1_sliding(:)) = 0;  solver.problem.ubw(ind_lambda1_sliding(:)) = 0;
solver.problem.lbw(ind_lambda0_sliding(:)) = 0;  solver.problem.ubw(ind_lambda0_sliding(:)) = 0;

%% Solve nlp
tic
[results,stats] = solver.solve();
cpu_time_iter = toc;
complementarity_iter = stats.complentarity_stats(end);
results = extract_results_from_solver(model,problem,settings,results);
results.complementarity_iter = complementarity_iter;
%% Verbose
if print_level>=2
    fprintf('-----------------------------------------------------------------------------------------------\n');
%     fprintf('Polishing step completed, complementarity resiudal %2.2e.\n',complementarity_iter);
    if model.n_u >0
        fprintf('CPU time of iteration: %2.2f s.\n',cpu_time_iter);
        fprintf('Objective function value: %2.4e.\n',cpu_time_iter);
        if time_optimal_problem
            fprintf('Final time T_opt: %2.4f.\n',w_opt(solver.problem.ind_t_final));
        end
    else
        fprintf('CPU time of iteration: %2.2f s.\n',cpu_time_iter);
    end
    fprintf('-----------------------------------------------------------------------------------------------\n');
end

end

