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

classdef Ipopt < handle % TODO maybe handle not necessary, revisit.
    properties

    end

    methods
        function solver = construct_solver(obj, nlp, solver_options, time_remaining)
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

            nlp.create_solver(opts_casadi_nlp);
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
            fprintf('\n|%-5s|%-10s|%-10s|%-10s|%-10s|%-10s|%-10s|%-10s|%-30s\n', 'iter', 'sigma', 'compl_res', 'inf_pr', 'inf_du', 'objective', 'CPU time', 'NLP iter', 'status');
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
            fprintf('|%-5d|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-10.3f|%-10d|%-30s\n',...
                ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                objective, stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status);
            
            % fprintf('%d\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.3f\t%d\t\t%s \n',...
            %     ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
            %     objective, stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status); 
        end
    end
end
