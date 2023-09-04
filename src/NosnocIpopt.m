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

classdef NosnocIpopt < handle % TODO maybe handle not necessary, revisit.
    properties

    end

    methods
        function solver = construct_solver(obj, nlp, solver_options, time_remaining)
            import casadi.*
            w = nlp.w;
            g = nlp.g;
            p = nlp.p;

            casadi_nlp = struct('f', nlp.augmented_objective, 'x', w, 'g', g, 'p', p);

            opts_casadi_nlp = solver_options.opts_casadi_nlp;
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'snopt');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'worhp');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'uno');
            if solver_options.timeout_cpu
                if exist('time_remaining')
                    opts_casadi_nlp.ipopt.max_cpu_time = time_remaining;
                else
                    opts_casadi_nlp.ipopt.max_cpu_time = solver_options.timeout_cpu;
                end
            elseif solver_options.timeout_wall
                if exist('time_remaining')
                    opts_casadi_nlp.ipopt.max_wall_time = time_remaining;
                else
                    opts_casadi_nlp.ipopt.max_wall_time = solver_options.timeout_wall;
                end
            end

            if ~solver_options.multiple_solvers
                solver = nlpsol(solver_options.solver_name, solver_options.solver, casadi_nlp, opts_casadi_nlp);
            else
                solver = {};
                sigma_k = solver_options.sigma_0;
                for k = 1:solver_options.N_homotopy
                    opts_casadi_nlp.ipopt.mu_init = sigma_k * 1e-1;
                    opts_casadi_nlp.ipopt.mu_target = sigma_k * 1e-1;
                    opts_casadi_nlp.ipopt.bound_relax_factor = sigma_k^2 * 1e-2;
                    opts_casadi_nlp.ipopt.mu_strategy = 'monotone';
                    if k == 1
                        opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
                        opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-4 * sigma_k;
                        opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-4 * sigma_k;
                    end
                    solver{k} = nlpsol(solver_options.solver_name, 'ipopt', casadi_nlp, opts_casadi_nlp);
                    % TODO: make homotopy update function and reuse here.
                    sigma_k = solver_options.homotopy_update_slope*sigma_k;
                end
            end
        end

        function solver_stats = cleanup_solver_stats(obj, solver_stats)
            if ~isfield(solver_stats, 'iterations')
                solver_stats.iterations = [];
            end
        end

        function failed = check_iteration_failed(obj, stats)
            switch stats.solver_stats(end).return_status
                case {'Solve_Succeeded', 'Solved_To_Acceptable_Level', 'Search_Direction_Becomes_Too_Small'}
                    failed = false;
                otherwise
                    failed = true;
            end
        end

        function timeout = check_timeout(obj, stats)
            switch stats.solver_stats(end).return_status
                case {'Maximum_WallTime_Exceeded', 'Maximum_CpuTime_Exceeded'}
                    timeout = 1;
                otherwise
                    timeout = 0;
            end
        end

        function w_opt = w_opt_from_results(obj, nlp_results)
            w_opt = full(nlp_results.x);
        end
        function f = f_from_results(obj, nlp_results)
            f = nlp_results.f;
        end        
        function g = g_from_results(obj, nlp_results)
            g = nlp_results.g;
        end

        function print_nlp_iter_header(obj)
            fprintf('\niter\t sigma \t\t compl_res\t inf_pr \t inf_du \t objective \t CPU time \t NLP iter\t status \n');
        end
        
        function print_nlp_iter_info(obj, stats)
            solver_stats = stats.solver_stats(end);
            ii = size(stats.solver_stats, 2);

            if isfield(solver_stats, 'iterations') && ~isempty(solver_stats.iterations)
                inf_pr = solver_stats.iterations.inf_pr(end);
                inf_du = solver_stats.iterations.inf_du(end);
                objective = solver_stats.iterations.obj(end);
            else
                inf_pr = nan;
                inf_du = nan;
                objective = nan;
            end
            fprintf('%d\t%6.2e\t %6.2e\t %6.2e\t %6.2e \t %6.2e \t %6.3f \t %d \t %s \n',...
                ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                objective, stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status);
        end
    end
end
