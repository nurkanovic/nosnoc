classdef NosnocModel < handle

    properties
        model_name = 'nosnoc_model'
        %----- basic user input -----
        % state
        x
        x0
        lbx
        ubx

        % user algebraics
        z
        z0
        lbz
        ubz
        g_z % uzer algebraic constraints

        % global variables
        v_global
        v0_global
        lbv_global
        ubv_global

        % control
        u
        u0
        lbu
        ubu

        % global parameters
        p_global
        p_global_val

        % time varying parameters
        p_time_var
        p_time_var_val

        % All params
        p

        F         % Dynamic functions
        c         % Switching functions
        S         % Sign matrix
        g_Stewart % Stewart indicator functions
        g_ind

        f_q        % Stage cost
        f_q_T      % Terminal cost

        % least squares
        lsq_x
        x_ref
        f_lsq_x
        x_ref_val
        lsq_u
        u_ref
        f_lsq_u
        u_ref_val
        lsq_T
        x_ref_end
        f_lsq_T
        x_ref_end_val

        % constraints
        g_path
        g_path_fun
        g_path_lb
        g_path_ub
        g_comp_path
        g_comp_path_fun
        g_terminal
        g_terminal_fun
        g_terminal_lb
        g_terminal_ub

        %----- DCS/time_freezing mode user input -----
        f_c  % Gap functions
        mu_f % Friction coef
        e    % Restitution coef
        q    % Generalized position
        v    % Generalized velocity
        f_v  % Generalized forces
        M    % Inertia matrix
        invM % Inertia matrix

        J_normal
        J_tangent
        D_tangent

        %----- Algebraic variables -----

        % all
        z_all

        % switching, only for step and stewart
        z_switching

        % impulse only for cls
        z_impulse

        % Stewart
        theta
        theta_sys
        lambda
        lambda_sys
        mu
        mu_sys

        % Step
        alpha
        alpha_sys
        lambda_n
        lambda_n_sys
        lambda_p
        lambda_p_sys
        beta
        gamma
        theta_step      % casadi expression either lifted or not
        theta_step_sys

        % CLS
        lambda_normal
        y_gap
        lambda_tangent
        gamma_d
        beta_d
        delta_d
        beta_conic
        gamma_conic
        p_vt
        n_vt
        alpha_vt
        % CLS impulse vars
        Lambda_normal
        Y_gap
        L_vn
        Lambda_tangent
        Gamma_d
        Beta_d
        Delta_d
        Gamma
        Beta
        P_vt
        N_vt
        Alpha_vt
        z_v

        %----- Generated model -----
        f_x % state dynamics (possibly not generated)

        g_switching
        g_lift
        g_z_all


        % Functions
        f_x_fun
        g_z_all_fun
        g_switching_fun
        c_fun
        dot_c_fun
        g_Stewart_fun
        g_impulse_fun
        f_c_fun
        M_fun
        invM_fun
        J_normal_fun
        J_tangent_fun
        D_tangent_fun
        f_q_T_fun
        J_cc_fun
        f_lsq_x_fun
        f_lsq_u_fun
        f_lsq_T_fun
        lambda00_fun
        g_lift_beta_fun
        g_lift_theta_step_fun

        % params
        p_time_var_stages
        p_dyn

        % flags
        friction_exists
        general_inclusion = 0
        there_exist_free_x0 = 0
        time_freezing_model_exists = 0
        custom_time_freezing = 0
        verified = 0
        vars_exist = 0
        equations_exist = 0

        % time freezing
        a_n
        k_aux

        % Dimensions
        dims
    end

    methods
        function obj = NosnocModel()
            check_matlab_requirement()
            obj.dims = NosnocDimensions();
        end

        function generate_equations(obj, problem_options)
            if obj.equations_exist
                return
            end
            import casadi.*
            dims = obj.dims;
            dcs_mode = problem_options.dcs_mode;

            %% Model functions of the DCS mode
            % if f_x doesnt exist we generate it from F
            % if it does we are in expert mode. TODO name.
            % Define differential equations
            if isempty(obj.f_x)
                obj.f_x = zeros(dims.n_x,1);
                % rhs of ODE;
                for ii = 1:dims.n_sys
                    switch dcs_mode
                      case 'Stewart'
                        obj.f_x = obj.f_x + obj.F{ii}*obj.theta_sys{ii};
                      case 'Heaviside'
                        obj.f_x = obj.f_x + obj.F{ii}*obj.theta_step_sys{ii};
                      case 'CLS'
                        if ~problem_options.lift_velocity_state
                            if obj.friction_exists
                                switch problem_options.friction_model
                                  case 'Conic'
                                    F_v = inv(obj.M)*(obj.f_v+obj.J_normal*obj.lambda_normal+obj.J_tangent*obj.lambda_tangent);
                                  case 'Polyhedral'
                                    F_v = inv(obj.M)*(obj.f_v+obj.J_normal*obj.lambda_normal+obj.D_tangent*obj.lambda_tangent);
                                end
                            else
                                F_v = inv(obj.M)*(obj.f_v+obj.J_normal*obj.lambda_normal);
                            end
                            obj.f_x = [obj.v;F_v];
                        else
                            obj.f_x = [v;obj.z_v];
                            if obj.friction_exists
                                switch problem_options.friction_model
                                  case 'Conic'
                                    g_lift_v = obj.M*obj.z_v -(obj.f_v +obj.J_normal*obj.lambda_normal + obj.J_tangent*obj.lambda_tangent);
                                  case 'Polyhedral'
                                    g_lift_v =  obj.M*obj.z_v -(obj.f_v + obj.J_normal*obj.lambda_normal + obj.D_tangent*obj.lambda_tangent);
                                end
                            else
                                g_lift_v =  obj.M*obj.z_v -(obj.f_v+obj.J_normal*obj.lambda_n);
                            end
                            obj.g_lift = [obj.g_lift;g_lift_v];
                        end
                    end
                end
            end

            %% Define algebraic equations

            g_switching = []; % collects switching function algebraic equations, 0 = g_i(x) - \lambda_i - e \mu_i, 0 = c(x)-lambda_p+lambda_n
            g_convex = []; % equation for the convex multiplers 1 = e' \theta
            g_alg_cls = []; % algebraic equations in a CLS
            f_comp_residual = 0; % the orthogonality cond)itions diag(\theta) \lambda = 0.
            g_impulse = [];
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                switch dcs_mode
                  case 'Stewart'
                    % basic algebraic equations and complementarity condtions of the DCS
                    % (Note that the cross complementarities are later defined when the discrete
                    % time variables for every IRK stage in the create_nlp_nosnoc function are defined.)
                    % g_ind_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_sys
                    % lambda_i'*theta_i = 0; for all i = 1,..., n_sys
                    % lambda_i >= 0;    for all i = 1,..., n_sys
                    % theta_i >= 0;     for all i = 1,..., n_sys
                    % Gradient of Lagrange Function of indicator LP
                    g_switching = [g_switching; obj.g_Stewart{ii}-obj.lambda_sys{ii}+obj.mu_sys{ii}*ones(dims.n_f_sys(ii),1)];
                    g_convex = [g_convex;ones(dims.n_f_sys(ii),1)'*obj.theta_sys{ii}-1];
                    lambda00_expr = [lambda00_expr; obj.g_Stewart{ii}- min(obj.g_Stewart{ii})];
                    f_comp_residual = f_comp_residual + obj.lambda_sys{ii}'*obj.theta_sys{ii};
                  case 'Heaviside'
                    % c_i(x) - (lambda_p_i-lambda_n_i)  = 0; for all i = 1,..., n_sys
                    % lambda_n_i'*alpha_i  = 0; for all i = 1,..., n_sys
                    % lambda_p_i'*(e-alpha_i)  = 0; for all i = 1,..., n_sys
                    % lambda_n_i >= 0;    for all i = 1,..., n_sys
                    % lambda_p_i >= 0;    for all i = 1,..., n_sys
                    % alpha_i >= 0;     for all i = 1,..., n_sys
                    g_switching = [g_switching;obj.c{ii}-obj.lambda_p_sys{ii}+obj.lambda_n_sys{ii}];
                    f_comp_residual = f_comp_residual + obj.lambda_n_sys{ii}'*obj.alpha_sys{ii}+obj.lambda_p_sys{ii}'*(ones(dims.n_c_sys(ii),1)-obj.alpha_sys{ii});
                    %             lambda00_expr = [lambda00_expr; max(c{ii},0); -min(c{ii}, 0)];
                    lambda00_expr = [lambda00_expr; -min(obj.c{ii}, 0); max(obj.c{ii},0)];
                  case 'CLS'
                    % dumy variables for impact quations:
                    v_post_impact = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'v_post_impact',dims.n_q);
                    v_pre_impact = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'v_pre_impact',dims.n_q);
                    g_alg_cls = [g_alg_cls; obj.y_gap - obj.f_c];
                    g_impulse = [g_impulse; obj.M*(v_post_impact-v_pre_impact)-obj.J_normal*obj.Lambda_normal]; % TODO: can this be relaxed? velocity junction
                    g_impulse = [g_impulse; obj.Y_gap-obj.f_c];
                    % add state jump for every contact
                    for ii = 1:dims.n_contacts
                        %                 g_impulse = [g_impulse; P_vn(ii)-N_vn(ii) - J_normal(:,ii)'*(v_post_impact+e(ii)*v_pre_impact)];
                        g_impulse = [g_impulse; obj.L_vn(ii) - obj.J_normal(:,ii)'*(v_post_impact+obj.e(ii)*v_pre_impact)];
                    end
                    if obj.friction_exists
                        switch problem_options.friction_model
                          case 'Conic'
                            g_impulse(1:dims.n_q) =  obj.M*(v_post_impact-v_pre_impact)-obj.J_normal*obj.Lambda_normal-obj.J_tangent*obj.Lambda_tangent;
                            % algebraic and friction equations
                            for ii = 1:dims.n_contacts
                                ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                                g_alg_cls  = [g_alg_cls;-obj.J_tangent(:,ind_temp)'*obj.v - 2*obj.gamma(ii)*obj.lambda_tangent(ind_temp);...
                                    obj.beta(ii)-((obj.mu_f(ii)*obj.lambda_normal(ii))^2-norm(obj.lambda_tangent(ind_temp))^2)];
                                g_impulse = [g_impulse;
                                    -obj.J_tangent(:,ind_temp)'*v_post_impact - 2*obj.Gamma(ii)*obj.Lambda_tangent(ind_temp);...
                                    obj.Beta - ((obj.mu_f(ii)*obj.Lambda_normal(ii))^2- norm(obj.Lambda_tangent(ind_temp))^2)];

                                if ~isequal(problem_options.conic_model_switch_handling,'Plain')
                                    % equality constraints for pos and neg parts of the tangetial velocity
                                    g_alg_cls  = [g_alg_cls;obj.J_tangent(:,ind_temp)'*obj.v-(obj.p_vt(ind_temp)-obj.n_vt(ind_temp))];
                                    g_impulse = [g_impulse;obj.J_tangent(:,ind_temp)'*v_post_impact - (obj.P_vt(ind_temp)-obj.N_vt(ind_temp))];
                                end
                            end
                          case 'Polyhedral'
                            g_impulse(1:dims.n_q) = obj.M*(v_post_impact-v_pre_impact)-obj.J_normal*obj.Lambda_normal-obj.D_tangent*obj.Lambda_tangent;
                            % impulse lifting equations
                            for ii = 1:dims.n_contacts
                                ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                                g_alg_cls  = [g_alg_cls;obj.beta_d(ii)-(obj.mu_f(ii)*obj.lambda_normal(ii) - sum(obj.lambda_tangent(ind_temp)));...
                                    obj.delta_d(ind_temp) - (obj.D_tangent(:,ind_temp)'*obj.v + obj.gamma_d(ii))];
                                g_impulse = [g_impulse;obj.Beta_d - (obj.mu_f(ii)*obj.Lambda_normal(ii)-sum(obj.Lambda_tangent(ind_temp)));...
                                    obj.Delta_d(ind_temp)- (obj.D_tangent(:,ind_temp)'*v_post_impact + obj.Gamma_d(ii))];
                            end
                        end
                    end
                    % some CLS functions
                    obj.M_fun = Function('M_fun', {obj.x}, {obj.M});
                    obj.invM_fun = Function('invM_fun', {obj.x}, {obj.invM});
                    obj.f_c_fun = Function('f_c_fun', {obj.x}, {obj.f_c});
                    obj.J_normal_fun = Function('J_normal_fun', {obj.x}, {obj.J_normal});
                    if obj.friction_exists
                        if isequal(problem_options.friction_model,'Conic')
                            obj.J_tangent_fun = Function('J_tangent_fun', {obj.x}, {obj.J_tangent});
                        else
                            obj.D_tangent_fun = Function('D_tangent_fun', {obj.x}, {obj.D_tangent});
                        end
                    end
                end
            end

            %% Lifting of forces in rigid bodies (either with time freezing or the dcs system)
            % the r.h.s of M(q)\ddot{q} = F(q,\dot{q},u) into  M{q}z_v-f(q,dor{q},u)=0; \ddot{q} = z_v
            % TODO: rename time_freezing_lift_forces into lift_velocity_state and
            % implemet above
            % g_lift_forces = [];
            % if time_freezing && time_freezing_lift_forces % TODO Is this broken with parameters/v_global
            %     f_v = f_x(n_q+1:2*n_q);
            %     if n_u > 0
            %         f_v_fun = Function('f_v_fun',{x,z_all,u,v_global},{f_v});
            %         z0_forces = full(f_v_fun(x0,z0_all,u0,v0_global));
            %     else % TODO remove this?
            %         f_v_fun = Function('f_v_fun',{x,z_all,v_global},{f_v});
            %         z0_forces = full(f_v_fun(x0,z0_all,v0_global));
            %     end
            %     z_forces = define_casadi_symbolic(casadi_symbolic_mode,'z_forces',n_q);
            %     z_all = [z_all;z_forces];
            %     n_z_all = n_z_all+n_q;
            %     z0_all = [z0_all;z0_forces];
            %     g_lift_forces = [M*z_forces - f_v]; % lifting function
            %     % new simple dynamics after lifting
            %     f_x = [f_x(1:n_q); z_forces; f_x(end-n_quad:end)];
            % end
            %%  collect all algebraic equations
            g_lp = [g_switching;g_convex;obj.g_lift];
            g_z_all = [g_lp; obj.g_z; g_alg_cls];
            n_algebraic_constraints = length(g_z_all);
            obj.g_switching = g_switching;
            obj.g_z_all = g_z_all;

            %% CasADi functions for indicator and region constraint functions
            % model equations
            % TODO: @ Anton make this to be a function of v_global as well (and test on telescop arm example)

            obj.f_x_fun = Function('f_x_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{obj.f_x,obj.f_q});
            obj.g_z_all_fun = Function('g_z_all_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{g_z_all}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
            obj.g_Stewart_fun = Function('g_Stewart_fun',{obj.x,obj.p},{obj.g_Stewart{:}});
            if isequal(dcs_mode,'CLS')
                obj.g_impulse_fun = Function('g_impulse_fun',{obj.q,v_post_impact,v_pre_impact,obj.z_impulse},{g_impulse});
            else
                c_all = vertcat(obj.c{:});
                obj.c_fun = Function('c_fun',{obj.x,obj.p},{c_all});
                dot_c = c_all.jacobian(obj.x)*obj.f_x;
                obj.dot_c_fun = Function('dot_c_fun',{obj.x,obj.z_all,obj.u,obj.p, obj.v_global},{dot_c}); % total time derivative of switching functions
                obj.lambda00_fun = Function('lambda00_fun',{obj.x,obj.p_global},{lambda00_expr});
                obj.g_switching_fun = Function('g_switching_fun', {obj.x,obj.z_switching,obj.u,obj.p}, {g_switching});
            end

            obj.J_cc_fun = Function('J_cc_fun',{obj.z_all},{f_comp_residual});
            obj.f_q_T_fun = Function('f_q_T',{obj.x,obj.p,obj.v_global},{obj.f_q_T});

            %%  CasADi functions for lest-square objective function terms
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{obj.x,obj.x_ref,obj.p},{obj.f_lsq_x});
            if dims.n_u > 0
                obj.f_lsq_u_fun = Function('f_lsq_u_fun',{obj.u,obj.u_ref,obj.p},{obj.f_lsq_u});
            end
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{obj.x,obj.x_ref_end,obj.p_global},{obj.f_lsq_T});

            if size(obj.g_terminal,1) ~= 0
                obj.g_terminal_fun  = Function('g_terminal_fun',{obj.x,obj.p_global,obj.v_global},{obj.g_terminal});
            end
            if size(obj.g_path,1) ~= 0
                obj.g_path_fun  = Function('g_path_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{obj.g_path});
            end
            if size(obj.g_comp_path,1) ~= 0
                obj.g_comp_path_fun  = Function('g_comp_path_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{obj.g_comp_path});
            end

            obj.equations_exist = 1;
        end

        function generate_variables(obj,problem_options)
            if obj.vars_exist
                return
            end
            import casadi.*
            casadi_symbolic_mode = problem_options.casadi_symbolic_mode;
            dcs_mode = problem_options.dcs_mode;
            dims = obj.dims;

            g_lift_theta_step = [];
            g_lift_beta = [];
            switch dcs_mode
              case 'Stewart'
                % dimensions
                dims.n_theta = sum(obj.dims.n_f_sys); % number of modes
                dims.n_lambda = dims.n_theta;
                for ii = 1:dims.n_sys
                    ii_str = num2str(ii);
                    % define theta (Filippov multiplers)
                    obj.theta_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['theta_' ii_str],obj.dims.n_f_sys(ii));
                    obj.theta = [obj.theta;obj.theta_sys{ii}];
                    % define mu_i (Lagrange multipler of e'theta =1;)
                    obj.mu_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['mu_' ii_str],1);
                    obj.mu = [obj.mu;obj.mu_sys{ii}];
                    % define lambda_i (Lagrange multipler of theta >= 0;)
                    obj.lambda_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_' ii_str],obj.dims.n_f_sys(ii));
                    obj.lambda = [obj.lambda;obj.lambda_sys{ii}];
                end
              case 'Heaviside'
                n_alpha = sum(obj.dims.n_c_sys);
                n_lambda_n = sum(obj.dims.n_c_sys);
                n_lambda_p = sum(obj.dims.n_c_sys);
                % for creae_nlp_fesd
                n_theta = 2*n_alpha;
                n_lambda = n_lambda_n+n_lambda_p;

                dims.n_alpha = n_alpha;
                dims.n_lambda_n = n_lambda_n;
                dims.n_lambda_p = n_lambda_p;
                dims.n_theta = n_theta;
                dims.n_lambda = n_lambda;

                % alpha
                if isempty(obj.alpha)
                    for ii=1:dims.n_sys
                        ii_str = num2str(ii);
                        % define alpha (selection of a set valued step function)
                        obj.alpha_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                        obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                    end
                else
                    idx = 1;
                    for ii=1:dims.n_sys
                        end_idx = idx + obj.dims.n_c_sys(ii)-1;
                        obj.alpha_sys{ii} = obj.alpha(idx:end_idx);
                        idx = end_idx+1;
                    end
                end
                

                if isempty(obj.lambda_n)
                    for ii=1:dims.n_sys
                        ii_str = num2str(ii);
                        % define lambda_n_i (Lagrange multipler of alpha >= 0;)
                        obj.lambda_n_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_n_' ii_str],obj.dims.n_c_sys(ii));
                        obj.lambda_n = [obj.lambda_n;obj.lambda_n_sys{ii}];
                    end
                else
                    idx = 1;
                    for ii=1:dims.n_sys
                        end_idx = idx + obj.dims.n_c_sys(ii)-1;
                        obj.lambda_n_sys{ii} = obj.lambda_n(idx:end_idx);
                        idx = end_idx+1;
                    end
                end

                if isempty(obj.lambda_p)
                    for ii=1:dims.n_sys
                        ii_str = num2str(ii);
                        % define lambda_p_i (Lagrange multipler of alpha <= 1;)
                        obj.lambda_p_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_p_' ii_str],obj.dims.n_c_sys(ii));
                        obj.lambda_p = [obj.lambda_p;obj.lambda_p_sys{ii}];
                    end
                else
                    idx = 1;
                    for ii=1:dims.n_sys
                        end_idx = idx + obj.dims.n_c_sys(ii)-1;
                        obj.lambda_p_sys{ii} = obj.lambda_p(idx:end_idx);
                        idx = end_idx+1;
                    end
                end

                if obj.general_inclusion % unpack alpha in this case
                    obj.alpha = vertcat(obj.alpha{:});
                end

                % Define already here lifting variables and functions
                % TODO allow for custom beta lifting
                % Theta collects the vector for dot_x = F(x)Theta,
                % terms or theta_step from lifting;
                if ~obj.general_inclusion
                    for ii = 1:dims.n_sys
                        theta_temp = [];
                        ii_str = num2str(ii);
                        S_temp = obj.S{ii};
                        if problem_options.pss_lift_step_functions
                            % TODO implement automatic lifting
                        else
                            if ~problem_options.time_freezing_inelastic
                                for j = 1:size(S_temp,1)
                                    alpha_ij = 1;
                                    for k = 1:size(S_temp,2)
                                        % create multiaffine term
                                        if S_temp(j,k) ~=0
                                            alpha_ij = alpha_ij * (0.5*(1-S_temp(j,k))+S_temp(j,k)*obj.alpha_sys{ii}(k) ) ;
                                        end
                                    end
                                    theta_temp = [theta_temp;alpha_ij];
                                end
                            end
                        end
                        obj.theta_step_sys{ii} = theta_temp;
                    end
                end

                %% time-freezing inelastic impacts (exploit structure with taiolored formulae)
                if problem_options.time_freezing_inelastic
                    % theta_step are the lifting variables that enter the ODE r.h.s.
                    if any(obj.mu_f > 0)
                        obj.friction_exists = 1;
                    else
                        obj.friction_exists = 0;
                      end

                    if ~problem_options.nonsmooth_switching_fun
                        alpha_q = obj.alpha(1:dims.n_contacts);
                        alpha_v_normal = obj.alpha(dims.n_contacts+1:2*dims.n_contacts);
                        if obj.friction_exists
                            alpha_v_tangent = obj.alpha(2*dims.n_contacts+1:end);
                        end
                    else
                        alpha_qv = obj.alpha(1:dims.n_contacts);
                        if obj.friction_exists
                            alpha_v_tangent = obj.alpha(dims.n_contacts+1:end);
                        end
                    end

                    obj.theta_step = define_casadi_symbolic(casadi_symbolic_mode,'theta_step',dims.n_aux+1);
                    theta_step_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_aux) '+1,1);']); % expressions for Filippov multipliers via alpha (and possibley beta).

                    % lifting variables
                    % empty expressions for initalization
                    beta_bilinear_ode = []; % for lifting bilinear terms in free flight dynamics multiplire
                    beta_bilinear_aux = []; % for lifiting bilinear terms appearing in aux. dynamics mutiplieres
                    beta_prod = []; % for lifting the multi affine term defineng the overall free flight dynamics multpliers
                                    % expressions for lifting
                    beta_bilinear_ode_expr = [];
                    beta_bilinear_aux_expr = [];
                    beta_prod_expr = [];
                    beta_prod_expr_guess = []; % extra expresion to make depend only on alpha (the one above depens on both and alpha and beta) - needed for eval. of inital guess

                    if problem_options.pss_lift_step_functions
                        % lift bilinear terms in product terms for free flight ode % (alpha_q*alpha_v)
                        if ~problem_options.nonsmooth_switching_fun
                            beta_bilinear_ode = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_ode',dims.n_contacts);
                            beta_bilinear_ode_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
                            if obj.friction_exists
                                % lift bilinear terms defining aux dynamics (1-alpha_q)*(1-alpha_v)
                                beta_bilinear_aux = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_aux',dims.n_contacts);
                                beta_bilinear_aux_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
                            end
                        end
                        if dims.n_contacts > 2
                            beta_prod = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear',dims.n_contacts-2);
                            beta_prod_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
                            beta_prod_expr_guess = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
                        end
                    end
                    obj.beta = [beta_bilinear_ode;
                        beta_bilinear_aux;
                        beta_prod];
                    % expresions for theta's and lifting
                    %% Filippov multipliers
                    alpha_ode = 1; % initalized product for free flight multiplier
                    if ~problem_options.pss_lift_step_functions
                        for ii = 1:dims.n_contacts
                            if problem_options.nonsmooth_switching_fun
                                alpha_ode = alpha_ode*alpha_qv(ii);
                                if obj.friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            else
                                alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
                                if obj.friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                                end
                            end
                        end
                        theta_step_expr(1) = alpha_ode;
                    else
                        % lift and have bilinear terms
                        if problem_options.nonsmooth_switching_fun
                            if dims.n_contacts <= 2
                                for ii = 1:dims.n_contacts
                                    alpha_ode = alpha_ode*alpha_qv(ii);
                                end
                                theta_step_expr(1) = alpha_ode;
                            else
                                beta_prod_expr(1) = (alpha_qv(1))*(alpha_qv(2));
                                beta_prod_expr_guess(1) = (alpha_qv(1))*(alpha_qv(2));
                                % lifting terms in between
                                for ii = 3:dims.n_contacts-1
                                    beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_qv(ii)); % beta_{i} = beta{i-1}*(prod_term_i+1}
                                    beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_qv(ii)); % this is to have an expression depending only on alpha for the inital guess eval
                                end
                                theta_step_expr(1)= beta_prod(end)*(alpha_qv(dims.n_contacts)); % last lifting term;
                            end
                            % lifting of aux dyn multiplier expressions
                            for ii = 1:dims.n_contacts
                                if obj.friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            end
                        else
                            % with two smooth switching functions
                            beta_bilinear_ode_expr = alpha_q.*alpha_v_normal;
                            if obj.friction_exists
                                beta_bilinear_aux_expr = (1-alpha_q).*(1-alpha_v_normal);
                            end

                            if dims.n_contacts <= 2
                                % here no lifting of product terms
                                for ii = 1:dims.n_contacts
                                    alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                end
                                theta_step_expr(1) = alpha_ode;
                            else
                                % here lifting of product terms
                                g_z_tf_beta_prod  = [beta_prod(1) - (alpha_q(1)+alpha_v_normal(1)-beta_bilinear_ode(1))*(alpha_q(2)+alpha_v_normal(2)-beta_bilinear_ode(2))]; % first lifting terms
                                                                                                                                                                              % lifting terms in between
                                for ii = 3:dims.n_contacts-1
                                    % beta_{i} = beta{i-1}*(prod_term_i+1}
                                    beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                    beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                end
                                % last lifting term;
                                theta_step_expr(1) = beta_prod(end)*(alpha_q(dims.n_contacts)+alpha_v_normal(dims.n_contacts)-beta_bilinear_ode(dims.n_contacts));
                            end
                            % lifting of aux dyn multiplier expressions
                            for ii = 1:dims.n_contacts
                                if obj.friction_exists
                                    theta_step_expr(ii+1) =  beta_bilinear_aux(ii)*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = beta_bilinear_aux(ii)*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                                end
                            end
                        end
                    end
                    % equality constraints in DCS
                    g_lift_theta_step = obj.theta_step-theta_step_expr;
                    g_lift_beta = obj.beta - [beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr];
                    % auxiliary functions to get inital guess for new algebraic variables theta and beta
                    obj.g_lift_theta_step_fun  = Function('g_lift_theta_step_fun',{obj.alpha,obj.beta},{theta_step_expr});
                    obj.g_lift_beta_fun = Function('g_lift_beta_fun',{obj.alpha},{[beta_bilinear_ode_expr;beta_bilinear_aux_expr;beta_prod_expr_guess]});
                    obj.theta_step_sys{1} = obj.theta_step;
                end
                dims.n_beta = length(obj.beta);
                dims.n_theta_step = length(obj.theta_step);
              case 'CLS'
                obj.lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'lambda_normal',dims.n_contacts);
                obj.y_gap = define_casadi_symbolic(casadi_symbolic_mode,'y_gap',dims.n_contacts);
                % Variables for impulse equations
                obj.Lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_normal',dims.n_contacts);
                obj.Y_gap = define_casadi_symbolic(casadi_symbolic_mode,'Y_gap',dims.n_contacts);
                %         P_vn = define_casadi_symbolic(casadi_symbolic_mode,'P_vn',dims.n_contacts); % pos part of state jump law
                %         N_vn = define_casadi_symbolic(casadi_symbolic_mode,'N_vn',dims.n_contacts); % neg part of state jump law
                obj.L_vn = define_casadi_symbolic(casadi_symbolic_mode,'L_vn',dims.n_contacts); % lifting variable for state jump law
                if obj.friction_exists
                    % tangetial contact froce (firction force)
                    obj.lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'lambda_tangent',dims.n_tangents);
                    % Impulse varaibles
                    obj.Lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_tangent',dims.n_tangents);
                    if isequal(problem_options.friction_model,'Polyhedral')
                        obj.gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'gamma_d',dims.n_contacts);
                        obj.beta_d = define_casadi_symbolic(casadi_symbolic_mode,'beta_d',dims.n_contacts); % lift friction cone bound
                        obj.delta_d = define_casadi_symbolic(casadi_symbolic_mode,'delta_d',dims.n_tangents); % lift lagrangian
                                                                                                     % Impulse varaibles
                        obj.Gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'Gamma_d',dims.n_contacts);
                        obj.Beta_d = define_casadi_symbolic(casadi_symbolic_mode,'Beta_d',dims.n_contacts); % lift friction cone bound
                        obj.Delta_d = define_casadi_symbolic(casadi_symbolic_mode,'Delta_d',dims.n_tangents); % lift lagrangian
                    end
                    if isequal(problem_options.friction_model,'Conic')
                        obj.gamma = define_casadi_symbolic(casadi_symbolic_mode,'gamma',dims.n_contacts);
                        obj.beta = define_casadi_symbolic(casadi_symbolic_mode,'beta',dims.n_contacts);
                        % Impulse variables;
                        obj.Gamma = define_casadi_symbolic(casadi_symbolic_mode,'Gamma',dims.n_contacts);
                        obj.Beta = define_casadi_symbolic(casadi_symbolic_mode,'Beta',dims.n_contacts);
                        switch problem_options.conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            obj.p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',dims.n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            obj.n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',dims.n_tangents); % negative parts of tagnetial velocity (for switch detection)
                                                                                                   % Impulse
                            obj.P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',dims.n_tangents);
                            obj.N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',dims.n_tangents);
                          case 'Lp'
                            obj.p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',dims.n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            obj.n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',dims.n_tangents); % negative parts of tagnetial velocity (for switch detection)
                            obj.alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'alpha_vt ',dims.n_tangents); % step function of tangential velocities
                                                                                                            % impulse
                            obj.P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',dims.n_tangents);
                            obj.N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',dims.n_tangents);
                            obj.Alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'Alpha_vt ',dims.n_tangents);
                        end
                    end
                end
            end
            obj.g_lift = [g_lift_theta_step;g_lift_beta];

            %% Collect algebaric varaibles for the specific DCS mode, define initial guess and bounds
            % TODO: @Anton: Do the bounds and guess specification already while defining the varaibles?

            switch dcs_mode
              case 'Stewart'
                % symbolic variables z = [theta;lambda;mu_Stewart];
                obj.z_all = [obj.theta;obj.lambda;obj.mu];
                obj.z_switching = [obj.lambda;obj.mu];
              case 'Heaviside'
                obj.z_all = [obj.alpha;obj.lambda_n;obj.lambda_p;obj.beta;obj.theta_step];
                obj.z_switching = [obj.lambda_n;obj.lambda_p];
              case 'CLS'
                obj.z_all = [obj.lambda_normal;obj.y_gap];
                % Impulse
                %         obj.z_impulse = [Lambda_normal;Y_gap;P_vn;N_vn];
                obj.z_impulse = [obj.Lambda_normal;obj.Y_gap;obj.L_vn];
                if obj.friction_exists
                    % tangetial contact froce (firction force)
                    obj.z_all = [obj.z_all;obj.lambda_tangent];
                    % Impulse
                    obj.z_impulse = [obj.z_impulse;obj.Lambda_tangent]; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
                    % friction aux multipliers
                    if isequal(problem_options.friction_model,'Polyhedral')
                        % polyhedral friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma_d;obj.beta_d;obj.delta_d];
                        % Polyhedral friction - collect impulse variables
                        obj.z_impulse = [obj.z_impulse;obj.Gamma_d;obj.Beta_d;obj.Delta_d];
                    end
                    if isequal(problem_options.friction_model,'Conic')
                        % conic friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma;obj.beta];
                        % Conic impulse
                        obj.z_impulse = [obj.z_impulse;obj.Gamma;obj.Beta];
                        switch problem_options.conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt];
                            % Impulse
                            obj.z_impulse = [obj.z_impulse;obj.P_vt;obj.N_vt];
                          case 'Lp'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt;obj.alpha_vt];
                            % Impulse
                            obj.z_impulse = [obj.z_impulse;obj.P_vt;obj.N_vt;obj.Alpha_vt];
                        end
                    end
                end
                if problem_options.lift_velocity_state
                    obj.z_v = define_casadi_symbolic(casadi_symbolic_mode,['z_v'],dims.n_q);
                    obj.z_all = [obj.z_all;obj.z_v];
                end
            end
            %% Add user provided algebraic
            obj.z_all = vertcat(obj.z_all,obj.z);

            obj.vars_exist = 1;
        end

        function verify_and_backfill(obj, problem_options)
            import casadi.*
            if problem_options.time_freezing
                obj.time_freezing(problem_options);
            end
            dims = obj.dims;

            if size(obj.x, 1) ~= 0
                dims.n_x = length(obj.x);
                % check  lbx
                if size(obj.lbx, 1) ~= 0
                    if length(obj.lbx) ~= dims.n_x
                        error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    obj.lbx = -inf*ones(dims.n_x,1);
                end
                % check ubx
                if size(obj.ubx, 1) ~= 0
                    if length(obj.ubx) ~= dims.n_x
                        error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    obj.ubx = inf*ones(dims.n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end
            problem_options.casadi_symbolic_mode = ['casadi.' obj.x(1).type_name()];

            %% Check is u provided
            if size(obj.u, 1) ~= 0
                dims.n_u = length(obj.u);
                % check  lbu
                if size(obj.lbu, 1) ~= 0
                    if length(obj.lbu) ~= dims.n_u
                        error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    obj.lbu = -inf*ones(dims.n_u,1);
                end
                % check ubu
                if size(obj.ubu, 1) ~= 0
                    if length(obj.ubu) ~= dims.n_u
                        error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    obj.ubu = inf*ones(dims.n_u,1);
                end
                % check u0
                if size(obj.u0, 1) ~= 0
                    if length(obj.u0) ~= dims.n_u
                        error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    obj.u0 = 0*ones(dims.n_u,1);
                end
            else
                obj.u = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if problem_options.print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
            end
            %% Check if z is provided
            if size(obj.z, 1) ~= 0
                dims.n_z = length(obj.z);

                if size(obj.z0, 1) ~= 0
                    if length(obj.z0) ~= dims.n_z
                        error('nosnoc: The vector z0, for the initial guess of z has the wrong size.')
                    end
                else
                    obj.z0 = zeros(dims.n_z, 1);
                end

                if size(obj.lbz, 1) ~= 0
                    if length(obj.lbz) ~= dims.n_z
                        error('nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.lbz = -inf*ones(dims.n_z, 1);
                end

                if size(obj.ubz, 1) ~= 0
                    if length(obj.ubz) ~= dims.n_z
                        error('nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.ubz = inf*ones(dims.n_z, 1);
                end
            else
                dims.n_z = 0;
                obj.z0 = [];
                obj.lbz = [];
                obj.ubz = [];
                obj.z = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'',0);
            end
            %% Global vars (i.e., variables that do not change with time)
            if size(obj.v_global, 1) ~= 0
                n_v_global = length(obj.v_global);
                if size(obj.v0_global, 1) ~= 0
                    if length(obj.v0_global) ~= n_v_global
                        error('nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
                    end
                else
                    obj.v0_global = zeros(n_v_global, 1);
                end

                if size(obj.lbv_global, 1) ~= 0
                    if length(obj.lbv_global) ~= n_v_global
                        error('nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
                    end
                else
                    obj.lbv_global = -inf*ones(n_v_global, 1);
                end

                if size(obj.ubv_global, 1) ~= 0
                    if length(obj.ubv_global) ~= n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(n_v_global, 1);
                end
            else
                n_v_global = 0;
                obj.v_global = define_casadi_symbolic(problem_options.casadi_symbolic_mode, '', 0);
                obj.v0_global = [];
                obj.lbv_global = [];
                obj.ubv_global = [];
            end

            %% Parameters (time variable and that do not change with time)
            if size(obj.p_global, 1) ~= 0
                dims.n_p_global = size(obj.p_global,1);
                if size(obj.p_global_val, 1) ~= 0
                    if size(obj.p_global_val,1) ~= dims.n_p_global
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(dims.n_p_global,1);
                end
            else
                dims.n_p_global = 0;
                obj.p_global = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, problem_options.N_stages]
                        error('nosnoc: User provided p_time_var_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, problem_options.N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:problem_options.N_stages
                    var_full = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['p_time_var_' num2str(ii)], dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'', [0, problem_options.N_stages]);
                obj.p_time_var_val = double.empty(0,problem_options.N_stages);
            end

            obj.p = vertcat(obj.p_global,obj.p_time_var);

            %% Stage and terminal costs check
            if size(obj.f_q, 1) == 0
                obj.f_q = 0;
            end

            if size(obj.f_q_T, 1) ~= 0
                terminal_cost = 1;
            else
                obj.f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if size(obj.lsq_x, 1) ~= 0
                if length(obj.lsq_x)<3
                    error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_x{2},1)~=size(obj.lsq_x{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(obj.lsq_x{1},1)~=size(obj.lsq_x{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
                end

                n_x_ref_rows = size(obj.lsq_x{2},1);
                n_x_ref_cols = size(obj.lsq_x{2},2);
                if n_x_ref_cols == problem_options.N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    obj.lsq_x{2} = repmat(obj.lsq_x{2},1,problem_options.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,problem_options.N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropriate size.')
                end
                obj.x_ref_val = obj.lsq_x{2};
                obj.x_ref = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                obj.f_lsq_x = (obj.lsq_x{1}-obj.x_ref)'*obj.lsq_x{3}*(obj.lsq_x{1}-obj.x_ref);
            else
                obj.x_ref = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'x_ref',1);
                obj.f_lsq_x = 0;
                obj.x_ref_val = zeros(1,problem_options.N_stages);
            end

            % least square terms for control inputs
            if size(obj.lsq_u, 1) ~= 0
                if length(obj.lsq_u)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_u{2},1)~=size(obj.lsq_u{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
                end
                if size(obj.lsq_u{1},1)~=size(obj.lsq_u{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
                end
                n_u_ref_rows = size(obj.lsq_u{2},1);
                n_u_ref_cols = size(obj.lsq_u{2},2);
                if n_u_ref_cols == problem_options.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    obj.lsq_u{2} = repmat(obj.lsq_u{2},1,problem_options.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,problem_options.N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropriate size.')
                end
                obj.u_ref_val = obj.lsq_u{2};
                obj.u_ref = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                obj.f_lsq_u = (obj.lsq_u{1}-obj.u_ref)'*obj.lsq_u{3}*(obj.lsq_u{1}-obj.u_ref);
            else
                obj.u_ref = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'u_ref',1);
                obj.f_lsq_u = 0;
                obj.u_ref_val = zeros(1,problem_options.N_stages);
            end


            % least square terms for control inputs
            if size(obj.lsq_T, 1) ~= 0
                % sanity chkecs on the input
                if length(obj.lsq_T)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference or the weight matrix are missing.')
                end
                if size(obj.lsq_T{2},1)~=size(obj.lsq_T{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
                end
                if size(obj.lsq_T{1},1)~=size(obj.lsq_T{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the terminal cost do not match.')
                end

                n_x_T_rows = size(obj.lsq_T{2},1);
                n_x_T_cols = size(obj.lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(obj.lsq_T{1}));
                    error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropriate size.')
                end
                obj.x_ref_end_val = obj.lsq_T{2};
                obj.x_ref_end = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                obj.f_lsq_T = (obj.lsq_T{1}-obj.x_ref_end)'*obj.lsq_T{3}*(obj.lsq_T{1}-obj.x_ref_end);
            else
                obj.x_ref_end  = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'x_ref_end',1);
                obj.f_lsq_T = 0;
                obj.x_ref_end_val = 0;
            end

            %% Inequality constraints check
            if size(obj.g_path, 1) ~= 0
                g_path_constraint  = 1;
                n_g_path = length(obj.g_path);
                if size(obj.g_path_lb, 1) ~= 0
                    if length(obj.g_path_lb)~=n_g_path;
                        error('The user provided vector g_path_lb has the wrong size.')
                    end
                else
                    obj.g_path_lb = -inf*ones(n_g_path,1);
                end

                if size(obj.g_path_ub, 1) ~= 0
                    if length(obj.g_path_ub)~=n_g_path;
                        error('The user provided vector g_path_ub has the wrong size.')
                    end
                else
                    obj.g_path_ub =  0*ones(n_g_path,1);
                end
            else
                n_g_path = 0;
                g_path_constraint  = 0;
            end

            %% Check path complementarity constraints
            g_comp_path_constraint  = 0;
            if size(obj.g_comp_path, 1) ~= 0
                g_comp_path_constraint = 1;
                if size(obj.g_comp_path, 2) ~= 2
                    error('g_comp_path must be of size (m, 2)')
                end
            else
                g_comp_path_constraint = 0;
            end
            %% Terminal constraints
            if size(obj.g_terminal, 1) ~= 0
                n_g_terminal = length(obj.g_terminal);
                if size(obj.g_terminal_lb, 1) ~= 0
                    if length(obj.g_terminal_lb)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_lb has the wrong size.')
                    end
                else
                    obj.g_terminal_lb = 0*ones(n_g_terminal,1);
                end

                if size(obj.g_terminal_ub, 1) ~= 0
                    if length(obj.g_terminal_ub)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_ub has the wrong size.')
                    end
                else
                    obj.g_terminal_ub =  0*ones(n_g_terminal,1);
                end
            else
                n_g_terminal = 0;
            end

            obj.g_Stewart = {};
            c_all = [];
            obj.friction_exists = 0;

            if isequal(problem_options.dcs_mode,'CLS')
                % TODO: there is some repetition to the time_freezing check, this should be unified!!!!
                % Check existence of relevant functions
                dims.n_sys = 1; % always one subystem in CLS (only loops over n_contacts later)
                if isempty(obj.f_c)
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                dims.n_contacts = length(obj.f_c);

                % coefficient of friction checks
                if size(obj.mu_f, 1) ~= 0
                    if length(obj.mu_f) ~= 1 && length(obj.mu_f) ~= dims.n_contacts
                        error('The length of model.mu_f has to be one or match the length of model.f_c')
                    end
                    if length(obj.mu_f) == 1
                        obj.mu_f = obj.mu_f*ones(dims.n_contacts,1);
                    end

                    if any(obj.mu_f > 0)
                        obj.friction_exists = 1;
                    else
                        obj.friction_exists = 0;
                    end
                else
                    obj.mu_f = zeros(dims.n_contacts,1);
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end
                if any(obj.mu_f<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                % coefficent of restiution check
                if isempty(obj.e)
                    error('nosnoc: Please provide a coefficient of restitution via model.e')
                else
                    if length(obj.e) ~= 1 && length(obj.e) ~= dims.n_contacts
                        error('The length of model.e has to be one or match the length of model.f_c')
                    end
                    if length(obj.e) == 1
                        obj.e = obj.e*ones(dims.n_contacts,1);
                    end
                end
                if any(abs(1-obj.e)>1) || any(obj.e<0)
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % dimensions and state space split
                if mod(dims.n_x,2)
                    dims.n_q = (dims.n_x-1)/2;
                else
                    dims.n_q = dims.n_x/2;
                end
                if isempty(obj.q) && isempty(obj.v)
                    obj.q = obj.x(1:dims.n_q);
                    obj.v = obj.x(dims.n_q+1:2*dims.n_q);
                    fprintf('Automatically detected state consisting of generalized coordinates q:\n');
                    print_casadi_vector(obj.q);
                    fprintf('and velocities v:\n');
                    print_casadi_vector(obj.v);
                end

                if isempty(obj.f_v)
                    error('nosnoc: the function f_v (collecting all generalized forces), in M(q) = dv/dt =  f_v(q,v,u) + J_n\lambda_n +J_t\lambda_t ~ is not provided in model.');
                end

                % Check intertia matrix
                if isempty(obj.M)
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    obj.M = eye(dims.n_q);
                    obj.invM = inv(obj.M);
                else
                    obj.invM = inv(obj.M);
                end

                %  Normal Contact Jacobian
                if size(obj.J_normal, 1) ~= 0
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(obj.J_normal,1)~=dims.n_q && size(obj.J_normal,2)~=dims.n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',dims.n_q,dims.n_contacts);
                        error('nosnoc: J_normal has the wrong size.')
                    end
                    J_normal_exists = 1;
                else
                    obj.J_normal = obj.f_c.jacobian(obj.q)';
                    fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                    J_normal_exists = 1;
                end

                % Tangent Contact Jacobian
                if obj.friction_exists
                    if isequal(problem_options.friction_model,'Conic')
                        if size(obj.J_tangent, 1) ~= 0
                            if size(obj.J_tangent,1)~=dims.n_q
                                error('nosnoc: J_tangent has the wrong size.')
                            end
                        else
                            error('nosnoc: please provide the tangent Jacobian in model.J_tangent.')
                        end
                    end

                    if isequal(problem_options.friction_model,'Polyhedral')
                        if isempty(obj.D_tangent)
                            error('nosnoc: please provide the polyhedral tangent Jacobian in model.D_tangent, e.g., using the conic tangent Jacobian model.J_tangent: D_tangent = [J_tangent(q_0),-J_tangent(q_0)].')
                        end
                    end
                end
                % Dimension of tangents
                dims.n_t = 0;
                if obj.friction_exists
                    if isequal(problem_options.friction_model,'Polyhedral')
                        dims.n_t = size(obj.D_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    elseif isequal(problem_options.friction_model,'Conic')
                        dims.n_t = size(obj.J_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    end
                    dims.n_tangents = dims.n_t*dims.n_contacts; % number tangent forces for all multpliers
                else
                    dims.n_tangents = 0;
                end
            end

            if isequal(problem_options.dcs_mode,'Heaviside') || isequal(problem_options.dcs_mode,'Stewart')
                if isempty(obj.F)
                    % Don't need F
                    if ~obj.general_inclusion
                        error('nosnoc: Matrix F (or matrices F_i) with PSS modes not provided.');
                    else
                        % TODO Implement more subsystems.
                        dims.n_sys = 1;
                    end
                else
                    % check how many subsystems are present
                    if iscell(obj.F)
                        dims.n_sys = length(obj.F);
                    else
                        obj.F = {obj.F};
                        dims.n_sys = 1;
                    end
                end

                if isempty(obj.S)
                    % if we are using general inclusions we dont need S.
                    if ~obj.general_inclusion
                        % if the matrix S is not provided, maybe the g_ind are available
                        % directly?
                        if isequal(problem_options.dcs_mode,'Stewart')
                            if ~isempty(obj.g_ind)
                                if ~iscell(obj.g_ind)
                                    obj.g_ind = {obj.g_ind};
                                end

                                for ii = 1:dims.n_sys
                                    % discriminant functions
                                    obj.g_Stewart{ii} = obj.g_ind{ii};
                                    obj.c{ii} = zeros(1,problem_options.casadi_symbolic_mode);
                                    c_all = [c_all; zeros(1,problem_options.casadi_symbolic_mode)];
                                end
                            else
                                error(['nosnoc: Neither the sign matrix S nor the indicator functions g_ind for regions are provided. ' ...
                                        'Either provide the matrix S and the expression for c, or the expression for g_ind.']);
                            end
                        else
                            error(['nosnoc: The user uses problem_options.dcs_mode = ''Heaviside'', but the sign matrix S is not provided. Please provide the matrix S and the expressions for c(x) (definfing the region boundaries).']);
                        end
                    else
                        if isempty(obj.c)
                            error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                        else
                            if ~iscell(obj.c)
                                obj.c = {obj.c};
                            end
                            if length(obj.c) ~= dims.n_sys
                                error('nosnoc: Number of different expressions for c does not match number of subsystems.')
                            end
                            for ii = 1:dims.n_sys
                                c_all = [c_all; obj.c{ii}];
                                n_c{ii} = length(obj.c{ii});
                                dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
                            end

                        end
                    end
                else
                    % Check if all data is avilable and if dimensions match.
                    if ~iscell(obj.S)
                        obj.S = {obj.S};
                    end
                    if length(obj.S) ~= dims.n_sys
                        error('nosnoc: Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
                    end
                    % Check constraint function c
                    if isempty(obj.c)
                        error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                    else
                        if ~iscell(obj.c)
                            obj.c = {obj.c};
                        end
                        if length(obj.c) ~= dims.n_sys
                            error('nosnoc: Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
                        end
                    end

                    % check are the matrices dense
                    if isequal(problem_options.dcs_mode,'Stewart')
                        for ii = 1:dims.n_sys
                            if any(sum(abs(obj.S{ii}),2)<size(obj.S{ii},2))
                                if dims.n_sys == 1
                                    error('nosnoc: The matrix S is not dense. Either provide a dense matrix or use problem_options.mode = ''Heaviside''.');
                                else
                                    error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use problem_options.mode = ''Heaviside''.']);
                                end
                            end
                        end
                    end

                    for ii = 1:dims.n_sys
                        if size(obj.S{ii},2) ~= length(obj.c{ii})
                            error('nosnoc: The matrix S and vector c do not have compatible dimension.');
                        end

                        % discrimnant functions
                        switch problem_options.dcs_mode
                            case 'Stewart'
                                % Create Stewart's indicator functions g_ind_ii
                                obj.g_Stewart{ii} = -obj.S{ii}*obj.c{ii};
                            case 'Heaviside'
                                %eval(['c_' num2str(ii) '= c{ii};']);
                        end
                        % dimensions of c
                        c_all = [c_all; obj.c{ii}];
                        n_c{ii} = length(obj.c{ii});
                        dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
                    end

                end

                if isempty(dims.n_c_sys)
                    dims.n_c_sys = 0;
                end

                if max(dims.n_c_sys) < 2 && isequal(problem_options.dcs_mode,'Heaviside')
                    pss_lift_step_functions = 0;
                    if problem_options.print_level >=1
                        fprintf('nosnoc: problem_options.pss_lift_step_functions set to 0, as are step fucntion selections are already entering the ODE linearly.\n')
                    end
                end

                if ~obj.general_inclusion
                    dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);
                else
                    dims.n_f_sys = [size(obj.f_x,1)];
                end
            end

            % populate dims
            obj.dims.n_s = problem_options.n_s;
        end

        function time_freezing(obj,problem_options)
            import casadi.*
            dims = obj.dims;
            if obj.custom_time_freezing
                return;
            end
            if ~isempty(obj.F)
                fprintf('nosnoc: model.F provided, the automated model reformulation will be not performed. \n')
                time_freezing_model_exists = 1;
            else
                time_freezing_model_exists = 0;
            end

            %% Experimental options
            inv_M_once = 0;
            %% Auxiliary functions
            sigma = SX.sym('sigma',1);
            a = SX.sym('a',1);
            b = SX.sym('b',1);
            f_natural_residual = 0.5*(b+a+sqrt((b-a+sigma)^2));
            % f_natural_residual = max(a,b);
            max_smooth_fun = Function('max_smooth_fun',{a,b,sigma},{f_natural_residual});
            %% Check is the provided user data valid and complete
            if ~time_freezing_model_exists
                % Check is there a gap gunctions
                if isempty(obj.f_c)
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                dims.n_contacts = length(obj.f_c);
                % check dimensions of contacts
                if isempty(dims.n_dim_contact)
                    warning('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
                    dims.n_dim_contact = 2;
                end

                % coefficent of restiution
                if isempty(obj.e)
                    error('nosnoc:  Please provide a coefficient of restitution via model.e')
                end

                if abs(1-obj.e)>1 || obj.e<0
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % coefficient of friction
                if isempty(obj.mu_f)
                    obj.mu_f = 0;
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end

                if length(obj.mu_f(:)) ~= 1 && length(obj.mu_f(:)) ~= dims.n_contacts
                    error('nosnoc: Vector mu has to have length 1 or n_contacts.')
                end

                if any(obj.mu_f<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                if any(obj.mu_f)>0
                    obj.friction_exists = 1;
                    if length(obj.mu_f(:)) == 1
                        obj.mu_f = ones(dims.n_contacts,1)*obj.mu_f;
                    end
                else
                    obj.friction_exists = 0;
                end

                % dimensions and state space split
                casadi_symbolic_mode = obj.x(1).type_name();
                if mod(size(obj.x,1),2)
                    dims.n_x = size(obj.x,1);
                    dims.n_q = (dims.n_x-1)/2;
                else
                    dims.n_x = size(obj.x,1);
                    dims.n_q = dims.n_x/2;
                end

                if isempty(obj.q) && isempty(obj.v)
                    obj.q = obj.x(1:dims.n_q);
                    obj.v = obj.x(dims.n_q+1:2*dims.n_q);
                end

                % check model function
                if isempty(obj.f_v)
                    error('nosnoc: the function f_v (collecting all generalized forces), in dv/dt =  f_v(q,v,u) + ... is not provided in model.');
                end

                % Check intertia matrix
                if isempty(obj.M)
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    obj.M = eye(dims.n_q);
                    obj.invM = inv(obj.M);
                else
                    obj.invM = inv(obj.M);
                end


                %  Normal Contact Jacobian
                if ~isempty(obj.J_normal)
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(obj.J_normal,1)~=dims.n_q && size(obj.J_normal,2)~=dims.n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',dims.n_q, dims.n_contacts);
                        error('nosnoc: J_normal has the wrong size.')
                    end
                    J_normal_exists = 1;
                else
                    obj.J_normal = obj.f_c.jacobian(obj.q)';
                    fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                    J_normal_exists = 1;
                end

                if is_zero(obj.J_normal)
                    error('nosnoc: The normal vector should have at least one non-zero entry.')
                end

                % Tangent Contact Jacobian
                if obj.friction_exists
                    if ~isempty(obj.J_tangent)
                        J_tangent_exists = 1;
                    else
                        J_tangent_exists = 0;
                    end

                    if J_tangent_exists
                        if size(obj.J_tangent,1)~=dims.n_q && size(obj.J_tangent,2)~=dims.n_contacts*(dims.n_dim_contact-1)
                            fprintf('nosnoc: J_tangent should be %d x %d matrix.\n',dims.n_q,dims.n_contacts*(dims.n_dim_contact-1));
                            error('nosnoc: J_tangent has the wrong size.')
                        end
                        J_tangent_exists = 1;
                    else
                        error('nosnoc: tangent Jacobian obj.J_tangent not provided.\n');
                    end
                else
                    J_tangent_exists = 0;
                end

                % qudrature state
                dims.n_quad  = 0;
                if problem_options.time_freezing_quadrature_state
                    % define quadrature state
                    L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
                    if ~isempty(obj.lbx)
                        obj.lbx = [obj.lbx;-inf];
                    end
                    if ~isempty(obj.ubx)
                        obj.ubx = [obj.ubx;inf];
                    end
                    obj.x = [obj.x;L];
                    obj.x0 = [obj.x0;0];
                    obj.f = [obj.f;obj.f_q];
                    obj.f_q = 0;
                    if ~isempty(obj.f_q_T)
                        obj.f_q_T  = obj.f_q_T + L;
                    else
                        obj.f_q_T = L;
                    end
                    dims.n_quad = 1;
                end
                % Clock state and dimensions
                if ~mod(dims.n_x,2)
                    % uneven number of states = it is assumed that the clock state is defined.
                    t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
                    % update lower and upper bounds of lbx and ubx
                    if ~isempty(obj.lbx)
                        obj.lbx = [obj.lbx;-inf];
                    end
                    if ~isempty(obj.ubx)
                        obj.ubx = [obj.ubx;inf];
                    end
                    obj.x = [obj.x;t];
                    obj.x0 = [obj.x0;0];
                end

                % normal and tangential velocities
                eps_t = 1e-7;
                v_normal = obj.J_normal'*obj.v;
                if obj.friction_exists
                    if dims.n_dim_contact == 2
                        v_tangent = (obj.J_tangent'*obj.v)';
                    else
                        v_tangent = obj.J_tangent'*obj.v;
                        v_tangent = reshape(v_tangent,2,dims.n_contacts); % 2 x n_c , the columns are the tangential velocities of the contact points

                    end
                    v_tangent_norms = [];
                    for ii = 1:dims.n_contacts
                        v_tangent_norms = [v_tangent_norms;norm(v_tangent(:,ii))];
                    end
                else
                    v_tangent  = [];
                end

                % parameter for auxiliary dynamics
                if isempty(obj.a_n)
                    obj.a_n  = 100;
                end
                %% Time-freezing reformulation
                if obj.e == 0
                    % Basic problem_options
                    problem_options.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
                    problem_options.dcs_mode = 'Heaviside'; % time freezing inelastic works better step (very inefficient with stewart)
                    %% switching function
                    if problem_options.nonsmooth_switching_fun
                        obj.c = [max_smooth_fun(obj.f_c,v_normal,0);v_tangent];
                    else
                        if dims.n_dim_contact == 2
                            obj.c = [obj.f_c;v_normal;v_tangent'];
                        else
                            obj.c = [obj.f_c;v_normal;v_tangent_norms-eps_t];
                        end
                    end
                    %% unconstrained dynamics with clock state
                    inv_M = inv(obj.M);
                    f_ode = [obj.v;...
                        inv_M*obj.f_v;
                        1];

                    %% Auxiliary dynamics
                    % where to use invM, in every aux dyn or only at the end
                    if inv_M_once
                        inv_M_aux = eye(dims.n_q);
                        inv_M_ext = blkdiag(zeros(dims.n_q),inv_M,0);
                    else
                        inv_M_aux = inv_M;
                        inv_M_ext = eye(dims.n_x+1);
                    end
                    f_aux_pos = []; % matrix wit all aux tan dyn
                    f_aux_neg = [];
                    % time freezing dynamics
                    if problem_options.stabilizing_q_dynamics
                        f_q_dynamics = -problem_options.kappa_stabilizing_q_dynamics*obj.J_normal*diag(obj.f_c);
                    else
                        f_q_dynamics = zeros(dims.n_q,dims.n_contacts);
                    end
                    f_aux_normal = [f_q_dynamics;inv_M_aux*obj.J_normal*obj.a_n;zeros(1,dims.n_contacts)];

                    for ii = 1:dims.n_contacts
                        if obj.friction_exists && obj.mu_f(ii)>0
                            % auxiliary tangent;
                            if dims.n_dim_contact == 2
                                v_tangent_ii = obj.J_tangent(:,ii)'*obj.v;
                                f_aux_pos_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(obj.J_normal(:,ii)-obj.J_tangent(:,ii)*(obj.mu_f(ii)))*obj.a_n;0]; % for v>0
                                f_aux_neg_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(obj.J_normal(:,ii)+obj.J_tangent(:,ii)*(obj.mu_f(ii)))*obj.a_n;0]; % for v<0
                            else
                                v_tangent_ii = v_tangent(:,ii);
                                f_aux_pos_ii = [f_q_dynamics(:,ii);inv_M_aux*(obj.J_normal(:,ii)*obj.a_n-obj.J_tangent(:,ii*2-1:ii*2)*obj.mu_f(ii)*obj.a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                                f_aux_neg_ii = [f_q_dynamics(:,ii);inv_M_aux*(obj.J_normal(:,ii)*obj.a_n+obj.J_tangent(:,ii*2-1:ii*2)*obj.mu_f(ii)*obj.a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                            end
                            f_aux_pos = [f_aux_pos,f_aux_pos_ii];
                            f_aux_neg= [f_aux_neg,f_aux_neg_ii];
                        end
                    end
                    % f_aux_normal = inv_M_aux*J_normal*a_n;
                    % f_aux_tangent = inv_M_aux*J_tangent*mu(ii)*a_n;
                    if obj.friction_exists
                        f_aux = [f_aux_pos,f_aux_neg];
                    else
                        f_aux = f_aux_normal;
                    end
                    obj.F = [f_ode (inv_M_ext*f_aux)];
                    obj.S = ones(size(obj.F,2),length(obj.c)); % dummy value to pass error checks
                                                   % number of auxiliary dynamicsm modes
                    if obj.friction_exists
                        dims.n_aux = 2*dims.n_contacts;
                    else
                        dims.n_aux = dims.n_contacts;
                    end
                else
                    % elastic
                    dcs_mode = 'Heaviside';
                    if isempty(obj.k_aux)
                        obj.k_aux = 10;
                        if problem_options.print_level > 1
                            fprintf('nosnoc: Setting default value for k_aux = 10.\n')
                        end
                    end
                    temp1 = 2*abs(log(obj.e));
                    temp2 = obj.k_aux/(pi^2+log(obj.e)^2);
                    c_aux = temp1/sqrt(temp2);
                    K = [0 1;-obj.k_aux -c_aux];
                    N  = [obj.J_normal zeros(dims.n_q,1);...
                        zeros(dims.n_q,1) obj.invM*obj.J_normal];
                    f_aux_n1 = N*K*N'*[obj.q;obj.v];
                    f_aux_n1 = [f_aux_n1;zeros(dims.n_quad+1,1)];
                    f_ode = [obj.v;obj.invM*obj.f_v;1];
                    % updated with clock state
                    obj.F = [f_ode, f_aux_n1];
                    obj.S = [1; -1];
                    obj.c = obj.f_c;
                    dims.n_aux = 1;
                end

                %% Settings updates
                obj.time_freezing_model_exists = 1;
                obj.dims.n_dim_contact = 2;
            end
        end

        function json = jsonencode(obj, varargin)
            import casadi.*
            model_struct = struct(obj);
            names = fieldnames(model_struct);
            for ii=1:numel(names)
                name = names{ii};
                if startsWith(class(model_struct.(name)), "casadi.")
                    model_struct.(name) = model_struct.(name).serialize();
                end
            end
            json = jsonencode(model_struct);
        end
    end  % methods

    methods(Static)
        function model = from_struct(model_struct)
            model = NosnocModel();

        end
    end
end % NosnocModel
