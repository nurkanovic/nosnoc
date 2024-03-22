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
    properties
        % General
        use_fesd(1,1) logical = 1
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX'

        % In case of simulation problem
        T_sim
        N_sim
        h_sim

        h % Step size
        h_k % Finite element step size
        T % Terminal time % TODO: why is there also T_val, do we need both?

        % descritization
        N_stages(1,1) {mustBeInteger, mustBePositive} = 1;
        N_finite_elements {mustBeInteger, mustBePositive} = 2;

        % IRK and FESD Settings
        n_s(1,1) {mustBeInteger} = 2
        irk_scheme(1,1) IRKSchemes = IRKSchemes.RADAU_IIA
        irk_representation IrkRepresentation = IrkRepresentation.integral;

        cross_comp_mode(1,1) CrossCompMode = CrossCompMode.FE_STAGE
        gamma_h(1,1) double {mustBeReal} = 1
        dcs_mode DcsMode = DcsMode.Stewart

        % lift complementarities
        lift_complementarities(1,1) logical = 0
        experimental_supervertical_form(1,1) logical = 0
        lower_bound_comp_lift(1,1) logical = 0

        % TODO: Initialization is really a solver concern, but there is not really a good extraction approach yet.
        % Initialization - Stewart
        lp_initialization(1,1) logical = 0

        % Initialization - Step
        initial_alpha(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_beta(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_theta_step(1,1) double {mustBeReal, mustBeFinite} = 1

        % Step theta lifting
        pss_lift_step_functions(1,1) logical = 0
        n_depth_step_lifting(1,1) {mustBeInteger, mustBeGreaterThanOrEqual(n_depth_step_lifting, 2)} = 2

        %General NLP/OCP Settings
        g_path_at_fe(1,1) logical = 0 % evaluate nonlinear path constraint at every finte element boundary
        g_path_at_stg(1,1) logical = 0 % evaluate nonlinear path constraint at every stage

        x_box_at_fe(1,1) logical = 1 % evaluate box constraint for diff states at every finite element boundary point
        x_box_at_stg(1,1) logical = 1 % evaluate box constraint for diff states at every stage point. (is set to zero per default in differential irk mode, as it becomes a linear instead of box constraint)
        time_optimal_problem(1,1) = 0

        % Step equilibration
        rho_h(1,1) double {mustBePositive} = 1
        step_equilibration(1,1) StepEquilibrationMode = StepEquilibrationMode.heuristic_mean % heuristic_mean, l2_relaxed, l2_relaxed_scaled, direct, direct_homotopy, direct_homotopy_lift
        step_equilibration_sigma(1,1) double {mustBePositive} = 0.1 % slope at zero in rescaling the indicator function, nu_ki_rescaled = tanh(nu_ki/step_equilibration_sigma)

        % Multiple shooting type discretization
        equidistant_control_grid(1,1) logical = 1

        % Time-Setting & Time-Freezing
        time_freezing(1,1) logical = 0
        time_freezing_inelastic(1,1) logical = 0
        % for time optimal problems and equidistant control grids in physical time
        use_speed_of_time_variables(1,1) logical = 0
        local_speed_of_time_variable(1,1) logical = 0
        stagewise_clock_constraint(1,1) logical = 1
        impose_terminal_phyisical_time(1,1) logical = 1
        s_sot0(1,1) double {mustBePositive} = 1
        s_sot_max(1,1) double {mustBePositive} = 25
        s_sot_min(1,1) double {mustBePositive} = 1
        S_sot_nominal(1,1) double {mustBePositive} = 1
        rho_sot(1,1) double {mustBeReal, mustBeNonnegative} = 0

        T_final_max(1,1) double {mustBePositive} = 1e2
        T_final_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        time_freezing_reduced_model(1,1) logical = 0 % analytic reduction of lifter formulation, less algebraic variables (experimental)
        time_freezing_hysteresis(1,1) logical = 0 % do not do automatic time freezing generation for hysteresis, it is not supported yet.
        time_freezing_nonlinear_friction_cone(1,1) logical = 1 % 1 - use nonlienar friction cone, 0 - use polyhedral l_inf approximation.

        time_freezing_quadrature_state(1,1) logical = 0 % make a nonsmooth quadrature state to integrate only if physical time is running
        time_freezing_lift_forces(1,1) logical = 0 % replace \dot{v} = M(q)^{-1}f(q,v,u) by dot{v} = z,  M(q)z - f(q,v,u) = 0

        % exerimental expert options
        nonsmooth_switching_fun(1,1) logical = 0 % experimental: use c = max(c1,c2) insetad of c = [c1c2]
        % expert mode: stabilize auxiliary dynamics in \nabla f_c(q) direction
        stabilizing_q_dynamics(1,1) logical = 0
        kappa_stabilizing_q_dynamics(1,1) double {mustBePositive} = 1e-5
        % Verbose
        print_level = 3

        % Settings specific to CLS
        friction_model (1,1) FrictionModel = FrictionModel.Conic; % use nonlinear friction ('Conic') or polyhedral approximation ('Polyhedral');
        conic_model_switch_handling (1,1) ConicModelSwitchHandling = ConicModelSwitchHandling.Abs; % How to treat switch detection in v_t.
        kappa_friction_reg (1,1) double {mustBeReal, mustBeNonnegative} = 0; % reg. term in friction equations to avoid large multipliers if no contact happens.
        lift_velocity_state(1,1) logical = 0; % define auxliary algebraic vairable, dot{v} = z_v, to avoid symbolic inversion of the inertia matrix;
        eps_cls double = 1e-3 % constraint: enforce f_c at Euler step with h * eps_cls


        % Relaxation of terminal constraint
        relax_terminal_constraint(1,1) {mustBeInteger} = 0; %  0  - hard constraint, 1 - ell_1 , 2  - ell_2 , 3 - ell_inf TODO enum
        relax_terminal_constraint_from_above(1,1) logical = 0;
        rho_terminal(1,1) double {mustBePositive} = 1e2;
        relax_terminal_constraint_homotopy(1,1) logical = 0; % terminal penalty is governed by homotopy parameter

        % Relaxation of terminal (or stage for equidistant grids) numerical/phyisical time constraints
        relax_terminal_numerical_time(1,1) logical = 0; % instead of imposing \sum_h = T, add it as ell1_penalty term
        rho_terminal_numerical_time(1,1) double {mustBeNonnegative} = 1e2
        relax_terminal_numerical_time_homotopy (1,1) logical = 0; % us the homotopy parameter for the penalty
        relax_terminal_physical_time(1,1) logical = 0; % instead of imposing \sum_h = T, add it as ell1_penalty term
        rho_terminal_physical_time(1,1) double {mustBeNonnegative} = 1e2
        relax_terminal_physical_time_homotopy (1,1) logical = 0; % us the homotopy parameter for the penalty

        % Experimental:
        no_initial_impacts(1,1) logical = 0
        cost_integration(1,1) logical = 1

        % All MPCC parameters
        T_val(1,1) double {mustBePositive} = 1
        p_val

        % Butcher Tableu
        A_irk double
        B_irk double
        b_irk double
        C_irk double
        D_irk double
        c_irk double

        right_boundary_point_explicit(1,1) logical
    end

    properties(Dependent)
        time_rescaling
    end

    methods
        function obj = Options()
            check_matlab_requirement()
            obj.p_val = [obj.rho_sot, obj.rho_h, obj.rho_terminal, obj.T_val];
        end

        function [] = preprocess(obj)
            import casadi.*

            % time grid
            % TODO: merge T_sim and T?
            if ~isempty(obj.N_sim) && ~isempty(obj.T_sim)
                obj.T = obj.T_sim/obj.N_sim;
                obj.h_sim = obj.T_sim/(obj.N_sim*obj.N_stages*obj.N_finite_elements);
                if obj.print_level >= 2 && exist("h_sim")
                    fprintf('Info: N_sim is given, so the h_sim provided by the user is overwritten.\n')
                end
            elseif ~isempty(obj.N_sim) || ~isempty(obj.T_sim)
                error('Provide both N_sim and T_sim for the integration.')
            end

            if numel(obj.T) ~= 1 && ~obj.time_optimal_problem
                error('terminal numerical time T must be provided if time_optimal_problem is False.');
            elseif numel(obj.T) == 0 && ~obj.time_optimal_problem
                obj.T = 1;
            elseif numel(obj.T) ~= 1
                error('terminal time T must be a positive scalar.');
            end
            obj.h = obj.T/obj.N_stages;
            obj.h_k = obj.h./obj.N_finite_elements;

            % check irk scheme compatibility
            if ismember(obj.irk_scheme, IRKSchemes.differential_only)
                if obj.print_level >=1
                    fprintf(['Info: The user provided RK scheme: ' char(obj.irk_scheme) ' is only available in the differential representation.\n']);
                end
                obj.irk_representation = 'differential';
            end
            if obj.n_s < 1 || obj.n_s > 9
                error("n_s must be in [1, 9]");
            end


            %
            if obj.relax_terminal_constraint > 3 || obj.relax_terminal_constraint < 0
                error("relax_terminal_constraint must be in [0, 3]");
            end
            if obj.gamma_h < 0 || obj.gamma_h > 1
                error("gamma_h must be in [0, 1]");
            end

            % check impacts mode
            if obj.dcs_mode == DcsMode.CLS && obj.time_freezing == 1
                if obj.print_level >=1
                    fprintf(['nosnoc: User uses dcs_mode.CLS and time_freezing = true at the same time. Defaulting to time freezing as it is currently more competitive\n']);
                end
                obj.time_freezing = 1;
                obj.dcs_mode = DcsMode.Step;
            end

            if obj.time_freezing
                obj.use_speed_of_time_variables = 1;
                obj.local_speed_of_time_variable = 1;
            end

            % N finite elements shape \
            if isscalar(obj.N_finite_elements)
                obj.N_finite_elements = obj.N_finite_elements*ones(obj.N_stages,1);
            else
                if length(obj.N_finite_elements) == obj.N_stages
                    obj.N_finite_elements = obj.N_finite_elements;
                else
                    error('settings.N_finite_elements must be length 1 or N_stages');
                end
            end

            if ~obj.time_rescaling
                if obj.print_level >= 1 && obj.use_speed_of_time_variables
                    warning('nosnoc:NosnocOptions:erroneous_use_speed_of_time_variables', "use_speed_of_time_variables erroneously set to true even though we are not rescaling time, this likely means your settings are somehow faulty. Using use_speed_of_time_variables = false")
                end
                obj.use_speed_of_time_variables = 0;
            end

            if ~obj.use_speed_of_time_variables
                if obj.print_level >= 1 && obj.local_speed_of_time_variable
                    warning('nosnoc:NosnocOptions:erroneous_local_speed_of_time_variable',"local_speed_of_time_variable erroneously set to true even though we are not using speed of time variales, this likely means your settings are somehow faulty. Using local_speed_of_time_variable = false")
                end
                obj.local_speed_of_time_variable = 0;
            end

            % create Butcher tableau
            % TODO this should live somewhere else. (i.e. butcher tableu should not be in settings)
            switch obj.irk_representation
              case IrkRepresentation.integral
                [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.n_s, obj.irk_scheme);
                if tau_root(end) == 1
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.B_irk = B;
                obj.C_irk = C;
                obj.D_irk = D;
                % also get time steps
                [~, ~, c_irk] = generate_butcher_tableu(obj.n_s,obj.irk_scheme);
                obj.c_irk = c_irk;

              case {IrkRepresentation.differential, IrkRepresentation.differential_lift_x}
                [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(obj.n_s,obj.irk_scheme);
                if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.A_irk = A_irk;
                obj.b_irk = b_irk;
                obj.c_irk = c_irk;
            end
            obj.right_boundary_point_explicit = right_boundary_point_explicit;
        end

        function time_rescaling = get.time_rescaling(obj)
            time_rescaling = (obj.time_freezing && obj.impose_terminal_phyisical_time) || obj.time_optimal_problem;
        end
    end
end
