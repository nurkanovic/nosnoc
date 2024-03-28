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

classdef NosnocIntegrator < handle
    properties
        model           % Nosnoc Model
        problem_options % Nosnoc Problem Options
        solver_options  % Nosnoc Solver options
        mpcc            % Nosnoc mpcc
        solver          % Nosnoc Solver

        u_sim           % simulated controls

        initial_guess   % Initial Guess for solver
    end

    methods
        function obj = NosnocIntegrator(model, problem_options, solver_options, u_sim, initial_guess)
            import casadi.*
            obj.model = model;
            obj.problem_options = problem_options;
            obj.solver_options = solver_options;
            obj.initial_guess = initial_guess;

            % update problem
            problem_options.equidistant_control_grid = 0;

            % Generate required components for the integrator
            obj.mpcc = NosnocMPCC(problem_options, model);
            obj.solver = NosnocSolver(obj.mpcc, solver_options);
            
            dims = model.dims;
            
            % TODO: u_sim is for sure broken. Just remove it and only use parameters? Dont need to be more convenient than casadi.
            % check does the provided u_sim has correct dims
            if dims.n_u > 0
                if ~exist('u_sim','var')
                    warning('no control values (u_sim) provided, using zeros')
                    u_sim = zeros(n_u, N_sim);
                else
                    [n_rows, n_cols] = size(u_sim);
                    if n_rows ~= dims.n_u || n_cols~=N_sim
                        error('Control inputs u_sim has the wrong size. Required dimension is n_u x N_sim.')
                    end
                end
            end
            obj.u_sim = u_sim;
        end

        function [results, integrator_stats] = solve(obj)
            solver = obj.solver;
            model = obj.model;
            mpcc = obj.mpcc;
            dims = model.dims;
            problem_options = obj.problem_options;
            solver_options = obj.solver_options;

            %% Initialization
            results = struct();
            names = get_result_names_from_settings(problem_options);
            names = [names, {"h"}];
            for name=names
                results.(name) = [];
                results.extended.(name) = [];
            end
            x0 = model.x0;
            results.x = x0;
            results.extended.x = x0;
            results.s_sot = [];
            if problem_options.dcs_mode == DcsMode.CLS
                results.x_with_impulse = x0;
                results.t_with_impulse = 0;
            end
            % stats
            complementarity_stats  = [];
            homotopy_iteration_stats = [];
            time_per_iter = [];
            sim_step_solver_results = [];
            t_current = 0;
            converged = [];
            constraint_violations = [];

            if ~isempty(obj.initial_guess)
                % remove duplicate indices
                [~, w] = unique(obj.initial_guess.t_grid, 'stable');
                % setdiff(w, 1:length(initial_guess.t_grid))
                % keyboard
                obj.initial_guess.t_grid = obj.initial_guess.t_grid(w);
                obj.initial_guess.x_traj = obj.initial_guess.x_traj(w, :);
                obj.initial_guess.lambda_normal_traj = obj.initial_guess.lambda_normal_traj(w, :);
            end


            %% Main simulation loop
            for ii = 1:problem_options.N_sim

                %% set problem parameters, controls
                solver.set("x0", x0);
                if dims.n_u > 0
                    solver.set('u', {u_sim(:,ii)});
                end

                if ii > 1 && solver_options.use_previous_solution_as_initial_guess
                    % TODO make this possible via solver interface directly
                    % (@anton) what does this look like
                    solver.mpcc.w0(dims.n_x+1:end) = res.w(dims.n_x+1:end);
                end

                %% set initial guess
                solver.set('x', {x0})
                solver.set('x_left_bp', {x0})
                if exist('initial_guess', 'var')
                    t_guess = t_current + cumsum([0; problem_options.h_k * ones(problem_options.N_finite_elements, 1)]);
                    x_guess = interp1(initial_guess.t_grid, initial_guess.x_traj, t_guess,'makima');
                    lambda_normal_guess = interp1(initial_guess.t_grid, initial_guess.lambda_normal_traj, t_guess(2:end-1), 'makima');
                    %
                    x_init = cell(1, problem_options.N_finite_elements);
                    y_gap_init = cell(1, problem_options.N_finite_elements);
                    for j = 1:problem_options.N_finite_elements
                        x_init{j} = x_guess(j+1, :);
                        y_gap_init{j} = full(model.f_c_fun(x_init{j}));
                    end
                    solver.set('x', x_init);
                    solver.set('x_left_bp', x_init);
                    %
                    if isequal(problem_options.dcs_mode, 'CLS')
                        ind_v = dims.n_q + 1: 2*dims.n_q;
                        solver.set('lambda_normal', {0});
                        % Note : set this to zero
                        L_vn_init = {0*(x_init{end}(ind_v(1)) + model.e * x_init{end-1}(ind_v(1)))};
                        solver.set('L_vn', L_vn_init);
                        %
                        diff_v = diff(x_guess(:,ind_v(1)));
                        Lambda_normal_init = max(abs(diff_v));
                        % Note: this is a bit hacky..
                        if Lambda_normal_init < 2.0
                            Lambda_normal_init = 0.0;
                        else
                            %                 keyboard
                        end
                        % Lambda_normal_init = lambda_normal_guess;
                        solver.set('Lambda_normal', {Lambda_normal_init});
                        % disp(['init Lambda_normal', num2str(lambda_normal_guess)]);
                        solver.set('y_gap', y_gap_init);
                    end
                end

                %% solve
                [res, solver_stats] = solver.solve();

                %% handle failure -> try second initialization
                if solver_stats.converged == 0
                    disp(['integrator_fesd: did not converge in step ', num2str(ii), ' constraint violation: ', num2str(solver_stats.constraint_violation, '%.2e')])
                    % TODO re-implement
                    %solver.print_iterate(res.W(:,end))
                    if problem_options.dcs_mode == "CLS"
                        disp('provided initial guess in integrator step did not converge, trying anther inital guess.');
                        solver.set('Lambda_normal', {7});
                        solver.set('lambda_normal', {0});
                        solver.set('y_gap', {0});
                        solver.set('Y_gap', {0});
                        solver.set('L_vn', {0});
                        [res, solver_stats] = solver.solve();
                        if solver_stats.converged == 0
                            disp(['integrator_fesd: did not converge in step ', num2str(ii), 'constraint violation: ', num2str(solver_stats.constraint_violation, '%.2e')])
                        end
                    end
                elseif solver_options.print_level >=2
                    fprintf('Integration step %d / %d (%2.3f s / %2.3f s) converged in %2.3f s. \n',...
                        ii, problem_options.N_sim, t_current, problem_options.T_sim, solver_stats.cpu_time_total);
                end

                %% gather results
                if solver_options.store_integrator_step_results
                    sim_step_solver_results = [sim_step_solver_results,res];
                end
                time_per_iter = [time_per_iter; solver_stats.cpu_time_total];
                constraint_violations = [constraint_violations, solver_stats.constraint_violation];
                converged = [converged, solver_stats.converged];

                %% Store data
                % update results struct
                for name=names
                    if name == 'x' || ~isfield(res, name)
                        continue
                    end
                    results.(name) = [results.(name), res.(name)];
                    if ~strcmp(name, 'h')
                        results.extended.(name) = [results.extended.(name), res.extended.(name)];
                    end
                end
                results.x = [results.x, res.x(:, 2:end)];
                results.extended.x = [results.extended.x, res.extended.x(:, 2:end)];

                % TODO: is there a better way to do this
                results.s_sot = [results.s_sot, res.w(flatten_ind(obj.mpcc.ind_sot))];

                if problem_options.dcs_mode == DcsMode.CLS
                    results.x_with_impulse = [results.x_with_impulse, res.x_with_impulse(:,2:end)];
                    % TODO maybe make solver take t0 and T_final as param?
                    results.t_with_impulse = [results.t_with_impulse, t_current + res.t_with_impulse(2:end)];
                end

                if solver_stats.converged == 0 && solver_options.break_simulation_if_infeasible
                    disp('solver did not converge and break_simulation_if_infeasible is True -> finish simulation with Failure!')
                    break
                end
                %% plot during execution
                if solver_options.real_time_plot
                    figure(100)
                    clf
                    grid on
                    hold on
                    xlim([0 T_sim]);
                    ylabel('$x(t)$','Interpreter','latex');
                    if problem_options.time_freezing
                        plot(results.x(end,:), results.x(1:end-1, :));
                        xlabel('$t$ [phyisical time]', 'Interpreter', 'latex');
                    else
                        t_temp = [0, cumsum(results.h)];
                        plot(t_temp, results.x(:, 1:end));
                        xlabel('$t$', 'Interpreter', 'latex');
                    end
                    ylim([min(min(results.x(1:end-1, :)))-0.3 max(max(results.x(1:end-1, :)))+0.3])
                end

                %% debug here before problem changes
                % if ~isequal(sign(x0(3:4)), sign(res.x(3:4,end)))
                %     disp('sign switch in velocity')
                %     keyboard
                % end

                %% update inital value
                x0 = res.x(:,end);
                t_current = t_current + problem_options.T;
                % update clock state
                if problem_options.impose_terminal_phyisical_time
                    solver.mpcc.p0(end) = solver.mpcc.p0(end)+problem_options.T;
                end
            end

            total_time = sum(time_per_iter);
            %% Verbose
            fprintf('\n');
            fprintf('-----------------------------------------------------------------\n');
            if problem_options.use_fesd
                fprintf( ['Simulation with the FESD ' char(problem_options.irk_scheme) ' with %d-RK stages completed.\n'], problem_options.n_s);
            else
                fprintf( ['Simulation with the standard ' char(problem_options.irk_scheme) ' with %d-RK stages completed.\n'], problem_options.n_s);
            end
            fprintf('---------------- Stats summary ----------------------------\n');
            fprintf('N_sim\t step-size\t\tN_stg\tN_FE\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter\tMax. comp.\tMin. comp.\n');
            fprintf('%d\t\t\t%2.3f\t\t%d\t\t%d\t\t%2.3f\t\t\t\t%2.3f\t\t\t%2.3f\t\t\t\t%2.2e\t%2.2e\n', problem_options.N_sim, problem_options.h_sim, problem_options.N_stages, problem_options.N_finite_elements(1), total_time, max(time_per_iter), min(time_per_iter), max(complementarity_stats), min(complementarity_stats));
            fprintf('-----------------------------------------------------------------\n\n');

            %% Output
            integrator_stats.complementarity_stats = complementarity_stats;
            integrator_stats.time_per_iter = time_per_iter;
            integrator_stats.homotopy_iteration_stats = homotopy_iteration_stats;
            integrator_stats.converged = converged;

            results.t_grid = cumsum([0,results.h])';

            % generate fine t_grid
            tgrid_long = 0;
            for ii = 1:length(results.h)
                for jj = 1:problem_options.n_s
                    tgrid_long = [tgrid_long; results.t_grid(ii) + problem_options.c_irk(jj)*results.h(ii)];
                end
            end
            results.extended.t_grid = tgrid_long;

            if solver_options.store_integrator_step_results
                results.sim_step_solver_results = sim_step_solver_results;
            end

            % Clear empty fields
            names = fieldnames(results);
            for k=1:numel(names)
                if isempty(results.(names{k}))
                    results = rmfield(results, names{k});
                end
            end
        end
    end
end
