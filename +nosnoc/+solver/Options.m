% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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
classdef Options < handle
% TODO clean up much of the work here.
    properties
        % General
        solver_name {mustBeTextScalar} = 'nosnoc_solver'
        solver = 'ipopt'
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX'

        % MPCC and Homotopy Settings
        complementarity_tol(1,1) double {mustBeReal, mustBePositive} = 1e-9
        objective_scaling_direct(1,1) logical = 1
        sigma_0(1,1) double {mustBeReal, mustBeNonnegative} = 1
        sigma_N(1,1) double {mustBeReal, mustBePositive} = 1e-9
        homotopy_update_rule = 'linear' % 'linear' sigma_k = homotopy_update_slope*sigma_N
                                        % 'superlinear' - sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent))
                                        % TODO enum
        assume_lower_bounds(1,1) logical = false;
        lift_complementarities(1,1) logical = false;
        
        homotopy_update_slope(1,1) double {mustBeReal} = 0.1
        homotopy_update_exponent(1,1) double {mustBeReal, mustBePositive} = 1.5 % the exponent in the superlinear rule
        N_homotopy = 0 % 0 -> set automatically
        s_elastic_max(1,1) double {mustBeReal, mustBePositive} = 1e1
        s_elastic_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        s_elastic_0(1,1) double {mustBeReal, mustBePositive} = 1
        decreasing_s_elastic_upper_bound(1,1) logical = 0

        polishing_step(1,1) logical = 0 % heuristic for fixing active set, yet exerimental, not recommended to use.
        polishing_derivative_test(1,1) logical = 0 % check in sliding mode also the derivative of switching functions

        % Verbose
        print_level = 3
        print_details_if_infeasible = 0;
        pause_homotopy_solver_if_infeasible = 0;

        % IPOPT Settings
        opts_casadi_nlp

        % Integrator Specific
        % TODO: Maybe a 3rd options set specifically for integrators :)
        use_previous_solution_as_initial_guess(1,1) logical = 0
        simulation_problem(1,1) logical = 0
        real_time_plot(1,1) logical = 0
        break_simulation_if_infeasible(1,1) logical = 0

        % TODO: make proper multiple solver class.
        multiple_solvers(1,1) logical = 0

        % All NLP parameters
        p_val

        %
        homotopy_steering_strategy(1,1) HomotopySteeringStrategy = HomotopySteeringStrategy.DIRECT
        lower_bound_relaxation(1,1) logical = 0

        % Output options
        store_integrator_step_results(1,1) logical = 0
        ipopt_callback = [] % This should be a function handle that takes (model,problem,settings,ipopt_solver,results)

        timeout_cpu(1,1) {mustBeReal, mustBeNonnegative} = 0;
        timeout_wall(1,1) {mustBeReal, mustBeNonnegative} = 0;

        warm_start_duals(1,1) logical = false;

        % Experimental
        normalize_homotopy_update(1,1) logical = 1
        norm_function
        calculate_stationarity_type(1,1) logical = 0;
    end

    properties(Dependent)
        time_rescaling
    end

    methods
        function obj = Options()

            default_tol = 1e-12;

            obj.opts_casadi_nlp.ipopt.print_level = 0;
            obj.opts_casadi_nlp.print_time = 0;
            obj.opts_casadi_nlp.ipopt.sb = 'yes';
            obj.opts_casadi_nlp.verbose = false;
            obj.opts_casadi_nlp.ipopt.max_iter = 500;
            obj.opts_casadi_nlp.ipopt.bound_relax_factor = 0;
            %obj.opts_casadi_nlp.ipopt.bound_relax_factor = 1e-8;
            %obj.opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
            obj.opts_casadi_nlp.ipopt.tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
            obj.opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
            obj.opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
            obj.opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
            %obj.opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-9;
            %obj.opts_casadi_nlp.ipopt.warm_start_bound_frac = 1e-9;
            % obj.opts_casadi_nlp.ipopt.warm_start_slack_bound_push = 1e-9;
            % obj.opts_casadi_nlp.ipopt.warm_start_slack_bound_frac = 1e-6;
            % obj.opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-9;
            % obj.opts_casadi_nlp.ipopt.warm_start_mult_init_max = 1e20;
            obj.opts_casadi_nlp.ipopt.linear_solver = 'mumps';
            obj.opts_casadi_nlp.snopt = struct();
            
            obj.opts_casadi_nlp.worhp = struct();
            obj.opts_casadi_nlp.uno = struct();

            obj.p_val = [obj.sigma_0];
        end

        function [] = preprocess(obj)
            import casadi.*
            % automatically set up casadi opts
            if obj.print_level < 4
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=0;
                obj.opts_casadi_nlp.ipopt.sb= 'yes';
                obj.opts_casadi_nlp.snopt.Minor_print_level = 0;
                obj.opts_casadi_nlp.snopt.Major_print_level = 0;
                %obj.opts_casadi_nlp.snopt.Solution = 'no';
                %obj.opts_casadi_nlp.snopt.Suppress_options_listings = 'no';
                %obj.opts_casadi_nlp.snopt.Summary_file = 0;
                
                
                obj.opts_casadi_nlp.worhp.NLPprint = -1;
                obj.opts_casadi_nlp.uno.statistics_print_header_every_iterations = '10000';
                
            elseif obj.print_level == 4
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=1;
                obj.opts_casadi_nlp.ipopt.sb= 'no';
                obj.opts_casadi_nlp.snopt.Minor_print_level = 1;
                obj.opts_casadi_nlp.snopt.Major_print_level = 1;
                obj.opts_casadi_nlp.worhp.NLPprint = 1;
            else
                obj.opts_casadi_nlp.ipopt.print_level = 5;
                obj.opts_casadi_nlp.worhp.NLPprint = 4;
                obj.opts_casadi_nlp.snopt.Minor_print_level = 11;
                obj.opts_casadi_nlp.snopt.Major_print_level = 1;

            end

            if any([obj.homotopy_update_slope >= 1, obj.homotopy_update_rule <= 0.0])
                error('homotopy_update_slope must be in (0, 1)');
            end

            obj.sigma_N = min(0.1*obj.complementarity_tol, obj.sigma_N); 
            if obj.N_homotopy == 0
                obj.N_homotopy = ceil(abs(log(obj.sigma_N / obj.sigma_0) / log(obj.homotopy_update_slope)));
                if ~strcmp(obj.homotopy_update_rule, 'linear')
                    warning('computing N_homotopy automatically only supported for linear homotopy_update_rule');
                end
            end
        end
    end
end
