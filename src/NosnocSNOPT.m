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

classdef NosnocSNOPT < handle % TODO maybe handle not necessary, revisit.
    properties

    end

    methods
        function solver = construct_solver(obj, nlp, solver_options, time_remaining)
            import casadi.*
            w = nlp.w;
            g = nlp.g;
            p = nlp.p;

            casadi_nlp = struct('f', nlp.cost, 'x', w, 'g', g, 'p', p);

            % TODO: Possible issue raise to casadi: allow unknown fields in options passed
            opts_casadi_nlp = solver_options.opts_casadi_nlp;
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'ipopt');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'worhp');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'uno');
            if solver_options.timeout_wall
                if exist('time_remaining')
                    opts_casadi_nlp.snopt.Time_limit = time_remaining;
                    opts_casadi_nlp.snopt.Timing_level = 3;
                else
                    opts_casadi_nlp.snopt.Time_limit = solver_options.timeout_wall;
                    opts_casadi_nlp.snopt.Timing_level = 3;
                end
            end
            solver = nlpsol(solver_options.solver_name, solver_options.solver, casadi_nlp, opts_casadi_nlp);
        end

        function solver_stats = cleanup_solver_stats(obj, solver_stats)
            if ~isfield(solver_stats, 'iterations')
                solver_stats.iterations = [];
            end
        end

        function failed = check_iteration_failed(obj, stats)
            switch stats.solver_stats(end).return_status
                case {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}
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
            f = full(nlp_results.f);
        end
        function g = g_from_results(obj, nlp_results)
            g = full(nlp_results.g);
        end
    end
end
