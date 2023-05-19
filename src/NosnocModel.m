classdef NosnocModel < handle

    properties
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
        
        f_q        % Stage cost
        f_q_T      % Terminal cost

        h % Step size
        h_k % Finite element step size

        N_finite_elements % Number of finite elements

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
        g_terminal
        g_terminal_fun
        g_terminal_lb
        g_terminal_ub

        % Terminal time
        T

        % Sim
        T_sim
        N_sim
        h_sim

        %----- DCS/time_freezing mode user input -----
        f_c  % Gap functions
        mu_f % Friction coef
        e    % Restitution coef
        q    % Generalized position
        v    % Generalized velocity
        f_v  % Generalized forces
        M    % Inertia matrix

        J_normal
        J_tangent
        D_tangent

        %----- Algebraic variables -----

        % all
        z_all

        % switching, only for step and stewart
        z_switching
        
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
        Lambda_tangent
        Gamma_d
        Beta_d
        Delta_d
        Gamma_conic
        Beta_conic
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
        
        % params
        p_time_var_stages
        p_dyn

        % flags
        friction_exists
        
        % Dimensions
        dims
    end

    methods
        function obj = NosnocModel()
            obj.dims = NosnocDimensions();
        end

        function reformulate(obj, settings)
            
        end

        function generate_diffeq(obj, settings)
            import casadi.*
            dims = obj.dims;
            dcs_mode = settings.dcs_mode;
            
            %% Model functions of the DCS mode
            % if f_x doesnt exist we generate it from F
            % if it does we are in expert mode. TODO name.
            % Define differential equations
            if isempty(obj.f_x)
                f_x = zeros(dims.n_x,1);
                % rhs of ODE;
                for ii = 1:dims.n_sys
                    switch dcs_mode
                      case 'Stewart'
                        f_x = f_x + obj.F{ii}*obj.theta_sys{ii};
                      case 'Step'
                        f_x = f_x + obj.F{ii}*obj.theta_step_sys{ii};
                      case 'CLS'
                        if ~lift_velocity_state
                            if settings.friction_exists
                                switch settings.friction_model
                                  case 'Conic'
                                    F_v = inv(M)*(f_v+obj.J_normal*obj.lambda_normal+obj.J_tangent*obj.lambda_tangent);
                                  case 'Polyhedral'
                                    F_v = inv(M)*(f_v+obj.J_normal*obj.lambda_normal+obj.D_tangent*obj.lambda_tangent);
                                end
                            else
                                F_v = inv(M)*(obj.f_v+obj.J_normal*obj.lambda_normal);
                            end
                            f_x = [v;F_v];
                        else
                            f_x = [v;obj.z_v];
                            if friction_exists
                                switch friction_model
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
            f_comp_residual = 0; % the orthogonality conditions diag(\theta) \lambda = 0.
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
                  case 'Step'
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
                    v_post_impact = define_casadi_symbolic(casadi_symbolic_mode,'v_post_impact',dims.n_q);
                    v_pre_impact = define_casadi_symbolic(casadi_symbolic_mode,'v_pre_impact',dims.n_q);
                    g_alg_cls = [g_alg_cls; y_gap - f_c];
                    g_impulse = [g_impulse; M*(v_post_impact-v_pre_impact)-J_normal*Lambda_normal]; % TODO: can this be relaxed? velocity junction
                    g_impulse = [g_impulse; Y_gap-f_c];
                    % add state jump for every contact
                    for ii = 1:dims.n_contacts
                        %                 g_impulse = [g_impulse; P_vn(ii)-N_vn(ii) - J_normal(:,ii)'*(v_post_impact+e(ii)*v_pre_impact)];
                        g_impulse = [g_impulse; L_vn(ii) - J_normal(:,ii)'*(v_post_impact+e(ii)*v_pre_impact)];
                    end
                    if friction_exists
                        switch friction_model
                          case 'Conic'
                            g_impulse(1:dims.n_q) =  M*(v_post_impact-v_pre_impact)-J_normal*Lambda_normal-J_tangent*Lambda_tangent;
                            % algebraic and friction equations
                            for ii = 1:dims.n_contacts
                                ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                                g_alg_cls  = [g_alg_cls;-J_tangent(:,ind_temp)'*v - 2*gamma(ii)*lambda_tangent(ind_temp);...
                                    beta(ii)-((mu(ii)*lambda_normal(ii))^2-norm(lambda_tangent(ind_temp))^2)];
                                g_impulse = [g_impulse;
                                    -J_tangent(:,ind_temp)'*v_post_impact - 2*Gamma(ii)*Lambda_tangent(ind_temp);...
                                    Beta - ((mu(ii)*Lambda_normal(ii))^2- norm(Lambda_tangent(ind_temp))^2)];

                                if ~isequal(conic_model_switch_handling,'Plain')
                                    % equality constraints for pos and neg parts of the tangetial velocity
                                    g_alg_cls  = [g_alg_cls;J_tangent(:,ind_temp)'*v-(p_vt(ind_temp)-n_vt(ind_temp))];
                                    g_impulse = [g_impulse;J_tangent(:,ind_temp)'*v_post_impact - (P_vt(ind_temp)-N_vt(ind_temp))];
                                end
                            end
                          case 'Polyhedral'
                            g_impulse(1:dims.n_q) = M*(v_post_impact-v_pre_impact)-J_normal*Lambda_normal-D_tangent*Lambda_tangent;
                            % impulse lifting equations
                            for ii = 1:dims.n_contacts
                                ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                                g_alg_cls  = [g_alg_cls;beta_d(ii)-(mu(ii)*lambda_normal(ii) - sum(lambda_tangent(ind_temp)));...
                                    delta_d(ind_temp) - (D_tangent(:,ind_temp)'*v + gamma_d(ii))];
                                g_impulse = [g_impulse;Beta_d - (mu(ii)*Lambda_normal(ii)-sum(Lambda_tangent(ind_temp)));...
                                    Delta_d(ind_temp)- (D_tangent(:,ind_temp)'*v_post_impact + Gamma_d(ii))];
                            end
                        end
                    end
                    % some CLS functions
                    M_fun = Function('M_fun', {x}, {M});
                    invM_fun = Function('invM_fun', {x}, {invM});
                    f_c_fun = Function('f_c_fun', {x}, {f_c});
                    J_normal_fun = Function('J_normal_fun', {x}, {J_normal});
                    if friction_exists
                        if isequal(friction_model,'Conic')
                            J_tangent_fun = Function('J_tangent_fun', {x}, {J_tangent});
                        else
                            D_tangent_fun = Function('D_tangent_fun', {x}, {D_tangent});
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

            %% CasADi functions for indicator and region constraint functions
            % model equations
            % TODO: @ Anton make this to be a function of v_global as well (and test on telescop arm example)

            obj.f_x_fun = Function('f_x_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{f_x,obj.f_q});
            obj.g_z_all_fun = Function('g_z_all_fun',{obj.x,obj.z_all,obj.u,obj.p,obj.v_global},{g_z_all}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
            obj.g_Stewart_fun = Function('g_Stewart_fun',{obj.x,obj.p},{obj.g_Stewart{:}});
            if isequal(dcs_mode,'CLS')
                obj.g_impulse_fun = Function('g_impulse_fun',{obj.q,v_post_impact,v_pre_impact,obj.z_impulse},{g_impulse});
            else
                c_all = vertcat(obj.c{:});
                obj.c_fun = Function('c_fun',{obj.x,obj.p},{c_all});
                dot_c = c_all.jacobian(obj.x)*f_x;
                obj.dot_c_fun = Function('c_fun',{obj.x,obj.z_all,obj.u,obj.p},{dot_c}); % total time derivative of switching functions
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
        end
        
        function generate_vars(obj,settings)
            import casadi.*
            casadi_symbolic_mode = settings.casadi_symbolic_mode;
            dcs_mode = settings.dcs_mode;
            dims = obj.dims;
            g_lift_theta_step = [];
            g_lift_beta = [];
            switch dcs_mode
              case 'Stewart'
                % dimensions
                dims.n_theta = sum(obj.dims.n_c_sys); % number of modes
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
              case 'Step'
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
                
                for ii = 1:dims.n_sys
                    ii_str = num2str(ii);
                    % define alpha (selection of a set valued step function)
                    if ~settings.general_inclusion
                        obj.alpha_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                        obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                    else
                        % TODO this needs to change if subsystems.
                        obj.alpha_sys{ii} = obj.alpha{ii};
                    end
                    % define lambda_n_i (Lagrange multipler of alpha >= 0;)
                    obj.lambda_n_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_n_' ii_str],obj.dims.n_c_sys(ii));
                    obj.lambda_n = [obj.lambda_n;obj.lambda_n_sys{ii}];
                    % define lambda_p_i (Lagrange multipler of alpha <= 1;)
                    obj.lambda_p_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_p_' ii_str],obj.dims.n_c_sys(ii));
                    obj.lambda_p = [obj.lambda_p;obj.lambda_p_sys{ii}];
                end

                % Define already here lifting variables and functions
                % TODO allow for custom beta lifting
                % Theta collects the vector for dot_x = F(x)Theta,
                % terms or theta_step from lifting;
                if ~settings.general_inclusion
                    for ii = 1:dims.n_sys
                        theta_temp = [];
                        ii_str = num2str(ii);
                        S_temp = obj.S{ii};
                        if settings.pss_lift_step_functions
                            % TODO implement automatic lifting
                        else
                            if ~settings.time_freezing_inelastic
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
                if settings.time_freezing_inelastic
                    % theta_step are the lifting variables that enter the ODE r.h.s.
                    if ~settings.nonsmooth_switching_fun
                        alpha_q = alpha(1:dims.n_contacts);
                        alpha_v_normal = alpha(dims.n_contacts+1:2*dims.n_contacts);
                        if friction_exists
                            alpha_v_tangent = obj.alpha(2*dims.n_contacts+1:end);
                        end
                    else
                        alpha_qv = alpha(1:dims.n_contacts);
                        if friction_exists
                            alpha_v_tangent = obj.alpha(dims.n_contacts+1:end);
                        end
                    end

                    obj.theta_step = define_casadi_symbolic(casadi_symbolic_mode,'theta_step',n_aux+1);
                    theta_step_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_aux) '+1,1);']); % expressions for Filippov multipliers via alpha (and possibley beta).

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

                    if pss_lift_step_functions
                        % lift bilinear terms in product terms for free flight ode % (alpha_q*alpha_v)
                        if ~nonsmooth_switching_fun
                            beta_bilinear_ode = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_ode',dims.n_contacts);
                            beta_bilinear_ode_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
                            if friction_exists
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
                    if ~pss_lift_step_functions
                        for ii = 1:dims.n_contacts
                            if nonsmooth_switching_fun
                                alpha_ode = alpha_ode*alpha_qv(ii);
                                if friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            else
                                alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
                                if friction_exists
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
                        if nonsmooth_switching_fun
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
                                if friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            end
                        else
                            % with two smooth switching functions
                            beta_bilinear_ode_expr = alpha_q.*alpha_v_normal;
                            if friction_exists
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
                                if friction_exists
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
                    g_lift_beta = beta - [beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr];
                    % auxiliary functions to get inital guess for new algebraic variables theta and beta
                    g_lift_theta_step_fun  = Function('g_lift_theta_step_fun',{alpha,beta},{theta_step_expr});
                    g_lift_beta_fun = Function('g_lift_beta_fun',{alpha},{[beta_bilinear_ode_expr;beta_bilinear_aux_expr;beta_prod_expr_guess]});
                    obj.theta_step_sys{1} = theta_step;
                end
                n_beta = length(obj.beta);
                n_theta_step = length(obj.theta_step);
              case 'CLS'
                lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'lambda_normal',n_contacts);
                y_gap = define_casadi_symbolic(casadi_symbolic_mode,'y_gap',n_contacts);
                % Variables for impulse equations
                Lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_normal',n_contacts);
                Y_gap = define_casadi_symbolic(casadi_symbolic_mode,'Y_gap',n_contacts);
                %         P_vn = define_casadi_symbolic(casadi_symbolic_mode,'P_vn',n_contacts); % pos part of state jump law
                %         N_vn = define_casadi_symbolic(casadi_symbolic_mode,'N_vn',n_contacts); % neg part of state jump law
                L_vn = define_casadi_symbolic(casadi_symbolic_mode,'L_vn',n_contacts); % lifting variable for state jump law
                if friction_exists
                    % tangetial contact froce (firction force)
                    lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'lambda_tangent',n_tangents);
                    % Impulse varaibles
                    Lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_tangent',n_tangents);
                    if isequal(friction_model,'Polyhedral')
                        gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'gamma_d',n_contacts);
                        beta_d = define_casadi_symbolic(casadi_symbolic_mode,'beta_d',n_contacts); % lift friction cone bound
                        delta_d = define_casadi_symbolic(casadi_symbolic_mode,'delta_d',n_tangents); % lift lagrangian
                                                                                                     % Impulse varaibles
                        Gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'Gamma_d',n_contacts);
                        Beta_d = define_casadi_symbolic(casadi_symbolic_mode,'Beta_d',n_contacts); % lift friction cone bound
                        Delta_d = define_casadi_symbolic(casadi_symbolic_mode,'Delta_d',n_tangents); % lift lagrangian
                    end
                    if isequal(friction_model,'Conic')
                        gamma = define_casadi_symbolic(casadi_symbolic_mode,'gamma',n_contacts);
                        beta = define_casadi_symbolic(casadi_symbolic_mode,'beta',n_contacts);
                        % Impulse variables;
                        Gamma = define_casadi_symbolic(casadi_symbolic_mode,'Gamma',n_contacts);
                        Beta = define_casadi_symbolic(casadi_symbolic_mode,'Beta',n_contacts);
                        switch conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                                                                                                   % Impulse
                            P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents);
                            N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents);
                          case 'Lp'
                            p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                            alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'alpha_vt ',n_tangents); % step function of tangential velocities
                                                                                                            % impulse
                            P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents);
                            N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents);
                            Alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'Alpha_vt ',n_tangents);
                        end
                    end
                end
            end
            obj.g_lift = [g_lift_theta_step;g_lift_beta];

            %% Collect algebaric varaibles for the specific DCS mode, define initial guess and bounds
            % TODO: @Anton: Do the bounds and guess specification already while defining the varaibles?

            z_impulse = []; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
            switch dcs_mode
              case 'Stewart'
                % symbolic variables z = [theta;lambda;mu_Stewart];
                obj.z_all = [obj.theta;obj.lambda;obj.mu];
                obj.z_switching = [obj.lambda;obj.mu];
              case 'Step'
                obj.z_all = [obj.alpha;obj.lambda_n;obj.lambda_p;obj.beta;obj.theta_step];
                obj.z_switching = [obj.lambda_n;obj.lambda_p];
              case 'CLS'
                obj.z_all = [obj.lambda_normal;obj.y_gap];
                % Impulse
                %         z_impulse = [Lambda_normal;Y_gap;P_vn;N_vn];
                z_impulse = [obj.Lambda_normal;obj.Y_gap;L_vn];
                if settings.friction_exists
                    % tangetial contact froce (firction force)
                    obj.z_all = [obj.z_all;obj.lambda_tangent];
                    % Impulse
                    z_impulse = [z_impulse;obj.Lambda_tangent]; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
                    % friction aux multipliers
                    if isequal(friction_model,'Polyhedral')
                        % polyhedral friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma_d;obj.beta_d;obj.delta_d];
                        % Polyhedral friction - collect impulse variables
                        z_impulse = [z_impulse;obj.Gamma_d;obj.Beta_d;obj.Delta_d];
                    end
                    if isequal(friction_model,'Conic')
                        % conic friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma;obj.beta];
                        % Conic impulse
                        z_impulse = [z_impulse;obj.Gamma;obj.Beta];
                        switch conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt];
                            % Impulse
                            z_impulse = [z_impulse;obj.P_vt;obj.N_vt];
                          case 'Lp'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt;obj.alpha_vt];
                            % Impulse
                            z_impulse = [z_impulse;obj.P_vt;obj.N_vt;obj.Alpha_vt];
                        end
                    end
                end
                if settings.lift_velocity_state
                    obj.z_v = define_casadi_symbolic(casadi_symbolic_mode,['z_v'],dims.n_q);
                    obj.z_all = [obj.z_all;obj.z_v];
                end
            end
            %% Add user provided algebraic
            obj.z_all = vertcat(obj.z_all,obj.z);
        end
        
        function verify_and_backfill(obj, settings)
            import casadi.*
            if settings.time_freezing
                obj.time_freezing(settings);
            end
            dims = obj.dims;

            if isempty(obj.dims.N_stages)
                if settings.print_level >= 2
                    fprintf('Info: number of stages N_stages not provided, setting to defaul value N_stages = 1.\n')
                end
                obj.dims.N_stages = 1;
            end

            obj.h = obj.T/dims.N_stages;
            % nominal lengths of the finite elements for different control intevrals, every control interval might have a different number of finite elements.
            dims.N_finite_elements = dims.N_finite_elements(:); % make a column vector of the input
            if length(dims.N_finite_elements) > dims.N_stages
                dims.N_finite_elements = dims.N_finite_elements(1:N_stages);
                if print_level >=1
                    fprintf('nosnoc: Provided N_finite_elements had more entries then N_stages, the surplus of entries was removed. \n')
                end
            end
            if length(dims.N_finite_elements) == 1
                dims.N_finite_elements = dims.N_finite_elements*ones(dims.N_stages,1);
            elseif length(dims.N_finite_elements) > 1 && length(dims.N_finite_elements) < dims.N_stages
                dims.N_finite_elements = dims.N_finite_elements(:); % make sure it is a column vector
                dims.N_finite_elements = [dims.N_finite_elements;dims.N_finite_elements(end)*ones(dims.N_stages-length(dims.N_finite_elements),1)];
            end
            obj.h_k = obj.h./dims.N_finite_elements;
            
            if ~isempty(obj.N_sim) && ~isempty(obj.T_sim)
                obj.T = obj.T_sim/obj.N_sim;
                obj.h_sim = obj.T_sim/(obj.N_sim*dims.N_stages*dims.N_finite_elements);
                if settings.print_level >= 2 && exist("h_sim")
                    fprintf('Info: N_sim is given, so the h_sim provided by the user is overwritten.\n')
                end
            elseif ~isempty(obj.N_sim) || ~isempty(obj.T_sim)
                error('Provide both N_sim and T_sim for the integration.')
            end

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
                obj.u = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if settings.print_level >=1
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
                obj.z = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
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
                    if length(obj.ubv_global) ~= dims.n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(n_v_global, 1);
                end
            else
                n_v_global = 0;
                obj.v_global = define_casadi_symbolic(settings.casadi_symbolic_mode, '', 0);
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
                obj.p_global = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
                if settings.print_level >= 1
                    fprintf('nosnoc: No global parameters given. \n')
                end
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, dims.N_stages]
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, dims.N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:dims.N_stages
                    var_full = define_casadi_symbolic(settings.casadi_symbolic_mode, ['p_time_var_' num2str(ii)], dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(settings.casadi_symbolic_mode,'', [0, dims.N_stages]);
                obj.p_time_var_val = double.empty(0,dims.N_stages);
                if settings.print_level >= 1
                    fprintf('nosnoc: No time varying parameters given. \n')
                end
            end

            obj.p = vertcat(obj.p_global,obj.p_time_var);

            %% Stage and terminal costs check
            if ~size(obj.f_q, 1) ~= 0
                if settings.print_level >=1
                    fprintf('nosnoc: No stage cost is provided. \n')
                end
                obj.f_q = 0;
            end

            if size(obj.f_q_T, 1) ~= 0
                terminal_cost = 1;
            else
                if settings.print_level >=1
                    fprintf('nosnoc: No terminal cost is provided. \n')
                end
                obj.f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if size(obj.lsq_x, 1) ~= 0
                if length(lsq_x)<3
                    error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(lsq_x{2},1)~=size(lsq_x{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(lsq_x{1},1)~=size(lsq_x{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
                end

                n_x_ref_rows = size(lsq_x{2},1);
                n_x_ref_cols = size(lsq_x{2},2);
                if n_x_ref_cols == dims.N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    lsq_x{2} = repmat(lsq_x{2},1,dims.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,dims.N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropaite size.')
                end
                obj.x_ref_val = lsq_x{2};
                obj.x_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                obj.f_lsq_x = (lsq_x{1}-x_ref)'*lsq_x{3}*(lsq_x{1}-x_ref);
            else
                obj.x_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref',1);
                obj.f_lsq_x = 0;
                obj.x_ref_val = zeros(1,dims.N_stages);
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
                if n_u_ref_cols == dims.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    lsq_u{2} = repmat(lsq_u{2},1,dims.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,dims.N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropaite size.')
                end
                obj.u_ref_val = lsq_u{2};
                obj.u_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                obj.f_lsq_u = (lsq_u{1}-u_ref)'*lsq_u{3}*(lsq_u{1}-u_ref);
            else
                obj.u_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'u_ref',1);
                obj.f_lsq_u = 0;
                obj.u_ref_val = zeros(1,dims.N_stages);
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

                n_x_T_rows = size(lsq_T{2},1);
                n_x_T_cols = size(lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(lsq_T{1}));
                    error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropaite size.')
                end
                obj.x_ref_end_val = lsq_T{2};
                obj.x_ref_end = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                obj.f_lsq_T = (lsq_T{1}-x_ref_end)'*lsq_T{3}*(lsq_T{1}-x_ref_end);
            else
                obj.x_ref_end  = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref_end',1);
                obj.f_lsq_T = 0;
                obj.x_ref_end_val = 0;
            end

            %% Inequality constraints check
            if size(obj.g_path, 1) ~= 0
                g_path_constraint  = 1;
                n_g_path = length(g_path);
                if size(obj.g_path_lb, 1) ~= 0
                    if length(g_path_lb)~=n_g_path;
                        error('The user provided vector g_path_lb has the wrong size.')
                    end
                else
                    obj.g_path_lb = -inf*ones(n_g_path,1);
                end

                if size(obj.g_path_ub, 1) ~= 0
                    if length(g_path_ub)~=n_g_path;
                        error('The user provided vector g_path_ub has the wrong size.')
                    end
                else
                    obj.g_path_ub =  0*ones(n_g_path,1);
                end
                obj.g_path_fun  = Function('g_path_fun',{obj.x,obj.u,obj.p,obj.v_global},{obj.g_path});
            else
                n_g_path = 0;
                g_path_constraint  = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No path constraints are provided. \n')
                end
            end

            %% Check path complementarity constraints
            g_comp_path_constraint  = 0;
            if size(obj.g_comp_path, 1) ~= 0
                g_comp_path_constraint  = 1;
                if size(g_comp_path, 2) ~= 2
                    error('g_comp_path must be of size (m, 2)')
                end
                obj.g_comp_path_fun  = Function('g_comp_path_fun',{obj.x,obj.u,obj.p,obj.v_global},{obj.g_comp_path});
            else
                g_comp_path_constraint = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No path complementarity constraints are provided. \n')
                end
            end
            %% Terminal constraints
            if size(obj.g_terminal, 1) ~= 0
                n_g_terminal = length(obj.g_terminal);
                if size(obj.g_terminal_lb, 1) ~= 0
                    if length(g_terminal_lb)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_lb has the wrong size.')
                    end
                else
                    obj.g_terminal_lb = 0*ones(n_g_terminal,1);
                end

                if size(obj.g_terminal_ub, 1) ~= 0
                    if length(g_terminal_ub)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_ub has the wrong size.')
                    end
                else
                    obj.g_terminal_ub =  0*ones(n_g_terminal,1);
                end
                obj.g_terminal_fun  = Function('g_terminal_fun',{obj.x,obj.p_global,obj.v_global},{obj.g_terminal});
            else
                n_g_terminal = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No terminal constraints are provided. \n')
                end
            end

            obj.g_Stewart = {};
            c_all = [];
            obj.friction_exists = 0;

            if isequal(settings.dcs_mode,'CLS')
                % TODO: there is some repetition to the time_freezing check, this should be unified!!!!
                % Check existence of relevant functions
                dims.n_sys = 1; % always one subystem in CLS (only loops over n_contacts later)
                if isempty(obj.f_c)
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                dims.n_contacts = length(obj.f_c);

                % coefficient of friction checks
                if size(obj.mu, 1) ~= 0
                    if length(obj.mu) ~= 1 && length(obj.mu) ~= dims.n_contacts
                        error('The length of model.mu has to be one or match the length of model.f_c')
                    end
                    if length(obj.mu) == 1
                        mu = mu*ones(dims.n_contacts,1);
                        obj.mu = mu;
                    end

                    if any(obj.mu > 0)
                        obj.friction_exists = 1;
                    else
                        obj.friction_exists = 0;
                    end
                else
                    obj.mu = zeros(dims.n_contacts,1);
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end
                if any(obj.mu<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                % coefficent of restiution check
                if isempty(obj.e)
                    error('nosnoc:  Please provide a coefficient of restitution via model.e')
                else
                    if length(obj.e) ~= 1 && length(obj.e) ~= dims.n_contacts
                        error('The length of model.e has to be one or match the length of model.f_c')
                    end
                    if length(obj.e) == 1
                        e = e*ones(dims.n_contacts,1);
                        obj.e = e;
                    end
                end
                if any(abs(1-e)>1) || any(e<0)
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % dimensions and state space split
                settings.casadi_symbolic_mode = obj.x(1).type_name();
                if mod(n_x,2)
                    dims.n_q = (n_x-1)/2;
                else
                    dims.n_q = n_x/2;
                end
                if isempty(obj.q) && isempty(obj.v)
                    q = x(1:dims.n_q);
                    v = x(dims.n_q+1:2*dims.n_q);
                end

                if isempty(obj.f_v)
                    error('nosnoc: the function f_v (collecting all generalized forces), in M(q) = dv/dt =  f_v(q,v,u) + J_n\lambda_n +J_t\lambda_t ~ is not provided in model.');
                end

                % Check intertia matrix
                if isempty(obj.M)
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    M = eye(dims.n_q);
                    invM = inv(M);
                else
                    invM = inv(M);
                end

                %  Normal Contact Jacobian
                if size(obj.J_normal, 1) ~= 0
                    J_normal = obj.J_normal;
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(J_normal,1)~=dims.n_q && size(J_normal,2)~=dims.n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',dims.n_q,dims.n_contacts);
                        error('nosnoc: J_normal has the wrong size.')
                    end
                    J_normal_exists = 1;
                else
                    J_normal = f_c.jacobian(q)';
                    fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                    J_normal_exists = 1;
                end

                if is_zero(J_normal)
                    error('nosnoc: The normal vector should have at least one non-zero entry.')
                end

                % Tangent Contact Jacobian
                if obj.friction_exists
                    if isequal(settings.friction_model,'Conic')
                        if size(obj.J_tangent, 1) ~= 0
                            J_tangent = obj.J_tangent;
                            if size(J_tangent,1)~=dims.n_q
                                error('nosnoc: J_tangent has the wrong size.')
                            end
                        else
                            error('nosnoc: please provide the tangent Jacobian in model.J_tangent.')
                        end
                    end

                    if isequal(settings.friction_model,'Polyhedral')
                        if isempty(obj.D_tangent)
                            error('nosnoc: please provide the polyhedral tangent Jacobian in model.D_tangent, e.g., using the conic tangent Jacobian model.J_tangent: D_tangent = [J_tangent(q_0),-J_tangent(q_0)].')
                        end
                    end
                end
                % Dimension of tangents
                dims.n_t = 0;
                if obj.friction_exists
                    if isequal(friction_model,'Polyhedral')
                        dims.n_t = size(D_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    elseif isequal(friction_model,'Conic')
                        dims.n_t = size(J_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    end
                    dims.n_tangents = dims.n_t*dims.n_contacts; % number tangent forces for all multpliers
                else
                    dims.n_tangents = 0;
                end
            end

            if isequal(settings.dcs_mode,'Step') || isequal(settings.dcs_mode,'Stewart')
                if isempty(obj.F)
                    % Don't need F
                    if ~settings.general_inclusion
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
                    if ~settings.general_inclusion
                        % if the matrix S is not provided, maybe the g_ind are available
                        % directly?
                        if isequal(settings.dcs_mode,'Stewart')
                            if exist('g_ind')
                                if ~iscell(g_ind)
                                    g_ind = {g_ind};
                                end

                                for ii = 1:dims.n_sys
                                    % discriminant functions
                                    obj.g_Stewart{ii} = g_ind{ii};
                                    c_all = [c_all; zeros(1,settings.casadi_symbolic_mode)];
                                end
                            else
                                error(['nosnoc: Neither the sign matrix S nor the indicator functions g_ind for regions are provided. ' ...
                                        'Either provide the matrix S and the expression for c, or the expression for g_ind.']);
                            end
                        else
                            error(['nosnoc: The user uses settings.dcs_mode = ''Step'', but the sign matrix S is not provided. Please provide the matrix S and the expressions for c(x) (definfing the region boundaries).']);
                        end
                    else
                        if isempty(obj.c)
                            error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                        else
                            if ~iscell(obj,c)
                                obj.c = {obj.c};
                            end
                            if length(obj.c) ~= dims.n_sys
                                error('nosnoc: Number of different expressions for c does not match number of subsystems.')
                            end
                            for ii = 1:dims.n_sys
                                c_all = [c_all; c{ii}];
                                n_c{ii} = length(c{ii});
                                dims.n_c_sys  = [dims.n_c_sys;length(c{ii})];
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
                    if isequal(settings.dcs_mode,'Stewart')
                        for ii = 1:dims.n_sys
                            if any(sum(abs(obj.S{ii}),2)<size(obj.S{ii},2))
                                if dims.n_sys == 1
                                    error('nosnoc: The matrix S is not dense. Either provide a dense matrix or use settings.mode = ''Step''.');
                                else
                                    error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use settings.mode = ''Step''.']);
                                end
                            end
                        end
                    end

                    for ii = 1:dims.n_sys
                        if size(obj.S{ii},2) ~= length(obj.c{ii})
                            error('nosnoc: The matrix S and vector c do not have compatible dimension.');
                        end

                        % discrimnant functions
                        switch settings.dcs_mode
                            case 'Stewart'
                                % Create Stewart's indicator functions g_ind_ii
                                obj.g_Stewart{ii} = -obj.S{ii}*obj.c{ii};
                            case 'Step'
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

                if max(dims.n_c_sys) < 2 && isequal(settings.dcs_mode,'Step')
                    pss_lift_step_functions = 0;
                    if settings.print_level >=1
                        fprintf('nosnoc: settings.pss_lift_step_functions set to 0, as are step fucntion selections are already entering the ODE linearly.\n')
                    end
                end

                if ~settings.general_inclusion
                    dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);
                else
                    dims.n_f_sys = [size(f_x,1)];
                end
            end

            % populate functions that can already be generated
            obj.c_fun = Function('c_fun',{obj.x,obj.p},{c_all});
            obj.g_Stewart_fun = Function('g_Stewart_fun',{obj.x,obj.p},{vertcat(obj.g_Stewart{:})});

            % populate dims
            obj.dims.n_s = settings.n_s;
        end

        function time_freezing(obj,settings)
            if isfield(model,'F')
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
            %% Chek is the provided user data valid and complete
            if ~time_freezing_model_exists
                % Check is there a gap gunctions
                if ~isfield(model,'f_c')
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                n_contacts = length(f_c);
                % check dimensions of contacts
                if ~isfield(model,'n_dim_contact')
                    warning('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
                    n_dim_contact = 2;
                end

                % coefficent of restiution
                if ~isfield(model,'e')
                    error('nosnoc:  Please provide a coefficient of restitution via model.e')
                end

                if abs(1-e)>1 || e<0
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % coefficient of friction
                if ~isfield(model,'mu')
                    mu = 0;
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end

                if length(mu(:)) ~= 1 && length(mu(:)) ~= n_contacts
                    errro('nosnoc: Vector mu has to have length 1 or n_c.')
                end

                if any(mu<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                if any(mu)>0
                    friction_exists = 1;
                    if length(mu(:)) == 1
                        mu = ones(n_contacts,1)*mu;
                    end
                else
                    friction_exists = 0;
                end

                % dimensions and state space split
                casadi_symbolic_mode = model.x(1).type_name();
                if mod(size(x,1),2)
                    n_x = size(x,1);
                    n_q = (n_x-1)/2;
                else
                    n_x = size(x,1);
                    n_q = n_x/2;
                end

                if ~isfield(model,'q') && ~isfield(model,'v')
                    q = x(1:n_q);
                    v = x(n_q+1:2*n_q);
                end

                % check model function
                if ~isfield(model,'f_v')
                    error('nosnoc: the function f_v (collecting all generalized forces), in dv/dt =  f_v(q,v,u) + ... is not provided in model.');
                end

                % Check intertia matrix
                if ~isfield(model,'M')
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    M = eye(n_q);
                    invM = inv(M);
                else
                    invM = inv(M);
                end


                %  Normal Contact Jacobian
                if isfield(model,'J_normal')
                    J_normal = model.J_normal;
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(J_normal,1)~=n_q && size(J_normal,2)~=n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',n_q,n_contacts);
                        error('nosnoc: J_normal has the wrong size.')
                    end
                    J_normal_exists = 1;
                else
                    J_normal = f_c.jacobian(q)';
                    fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                    J_normal_exists = 1;
                end

                if is_zero(J_normal)
                    error('nosnoc: The normal vector should have at least one non-zero entry.')
                end

                % Tangent Contact Jacobian
                if friction_exists
                    if isfield(model,'J_tangent')
                        J_tangent = model.J_tangent;
                        J_tangent_exists = 1;
                    else
                        J_tangent_exists = 0;
                    end

                    if J_tangent_exists
                        if size(J_tangent,1)~=n_q && size(J_tangent,2)~=n_contacts*(n_dim_contact-1)
                            fprintf('nosnoc: J_tangent should be %d x %d matrix.\n',n_q,n_contacts*(n_dim_contact-1));
                            error('nosnoc: J_tangent has the wrong size.')
                        end
                        J_tangent_exists = 1;
                    else
                        error('nosnoc: tangent Jacobian model.J_tangent not provided.\n');
                    end
                else
                    J_tangent_exists = 0;
                end

                % qudrature state
                n_quad  = 0;
                if settings.time_freezing_quadrature_state
                    % define quadrature state
                    L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
                    if isfield(model,'lbx')
                        model.lbx = [model.lbx;-inf];
                    end
                    if isfield(model,'ubx')
                        model.ubx = [model.ubx;inf];
                    end
                    x = [x;L];
                    model.x = x;
                    model.x0 = [model.x0;0];
                    f = [f;f_q];
                    model.f = f;
                    model.f_q = 0;
                    if isfield(model,'f_q_T')
                        model.f_q_T  = model.f_q_T + L;
                    else
                        model.f_q_T  = L;
                    end
                    n_quad = 1;
                end
                % Clock state and dimensions
                if ~mod(n_x,2)
                    % uneven number of states = it is assumed that the clock state is defined.
                    t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
                    % update lower and upper bounds of lbx and ubx
                    if isfield(model,'lbx')
                        model.lbx = [model.lbx;-inf];
                    end
                    if isfield(model,'ubx')
                        model.ubx = [model.ubx;inf];
                    end
                    x = [x;t];
                    x0 = [x0;0];
                end

                % normal and tangential velocities
                eps_t = 1e-7;
                v_normal = J_normal'*v;
                if friction_exists
                    if n_dim_contact == 2
                        v_tangent = (J_tangent'*v)';
                    else
                        v_tangent = J_tangent'*v;
                        v_tangent = reshape(v_tangent,2,n_contacts); % 2 x n_c , the columns are the tangential velocities of the contact points

                    end
                    v_tangent_norms = [];
                    for ii = 1:n_contacts
                        v_tangent_norms = [v_tangent_norms;norm(v_tangent(:,ii))];
                    end
                else
                    v_tangent  = [];
                end

                % parameter for auxiliary dynamics
                if ~isfield(model,'a_n')
                    a_n  = 100;
                end
                %% Time-freezing reformulation
                if e == 0
                    % Basic settings
                    settings.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
                    settings.dcs_mode = 'Step'; % time freezing inelastic works better step (very inefficient with stewart)
                    %% switching function
                    if settings.nonsmooth_switching_fun
                        c = [max_smooth_fun(f_c,v_normal,0);v_tangent];
                    else
                        if n_dim_contact == 2
                            c = [f_c;v_normal;v_tangent'];
                        else
                            c = [f_c;v_normal;v_tangent_norms-eps_t];
                        end
                    end
                    %% unconstrained dynamics with clock state
                    inv_M = inv(M);
                    f_ode = [v;...
                        inv_M*f_v;
                        1];

                    %% Auxiliary dynamics
                    % where to use invM, in every aux dyn or only at the end
                    if inv_M_once
                        inv_M_aux = eye(n_q);
                        inv_M_ext = blkdiag(zeros(n_q),inv_M,0);
                    else
                        inv_M_aux = inv_M;
                        inv_M_ext = eye(n_x+1);
                    end
                    f_aux_pos = []; % matrix wit all aux tan dyn
                    f_aux_neg = [];
                    % time freezing dynamics
                    if settings.stabilizing_q_dynamics
                        f_q_dynamics = -settings.kappa_stabilizing_q_dynamics*J_normal*diag(f_c);
                    else
                        f_q_dynamics = zeros(n_q,n_contacts);
                    end
                    f_aux_normal = [f_q_dynamics;inv_M_aux*J_normal*a_n;zeros(1,n_contacts)];

                    for ii = 1:n_contacts
                        if friction_exists && mu(ii)>0
                            % auxiliary tangent;
                            if n_dim_contact == 2
                                v_tangent_ii = J_tangent(:,ii)'*v;
                                f_aux_pos_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(J_normal(:,ii)-J_tangent(:,ii)*(mu(ii)))*a_n;0]; % for v>0
                                f_aux_neg_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(J_normal(:,ii)+J_tangent(:,ii)*(mu(ii)))*a_n;0]; % for v<0
                            else
                                v_tangent_ii = v_tangent(:,ii);
                                f_aux_pos_ii = [f_q_dynamics(:,ii);inv_M_aux*(J_normal(:,ii)*a_n-J_tangent(:,ii*2-1:ii*2)*mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                                f_aux_neg_ii = [f_q_dynamics(:,ii);inv_M_aux*(J_normal(:,ii)*a_n+J_tangent(:,ii*2-1:ii*2)*mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                            end
                            f_aux_pos = [f_aux_pos,f_aux_pos_ii];
                            f_aux_neg= [f_aux_neg,f_aux_neg_ii];
                        end
                    end
                    % f_aux_normal = inv_M_aux*J_normal*a_n;
                    % f_aux_tangent = inv_M_aux*J_tangent*mu(ii)*a_n;
                    if friction_exists
                        f_aux = [f_aux_pos,f_aux_neg];
                    else
                        f_aux = f_aux_normal;
                    end
                    F = [f_ode (inv_M_ext*f_aux)];
                    S = ones(size(F,2),length(c)); % dummy value to pass error checks
                                                   % number of auxiliary dynamicsm modes
                    if friction_exists
                        n_aux = 2*n_contacts;
                    else
                        n_aux = n_contacts;
                    end
                else
                    % elastic
                    dcs_mode = 'Step';
                    if ~isfield(model,'k_aux')
                        k_aux = 10;
                        if settings.print_level > 1
                            fprintf('nosnoc: Setting default value for k_aux = 10.\n')
                        end
                    end
                    temp1 = 2*abs(log(e));
                    temp2 = k_aux/(pi^2+log(e)^2);
                    c_aux = temp1/sqrt(temp2);
                    K = [0 1;-k_aux -c_aux];
                    N  = [J_normal zeros(n_q,1);...
                        zeros(n_q,1) invM*J_normal];
                    f_aux_n1 = N*K*N'*[q;v];
                    f_aux_n1 = [f_aux_n1;zeros(n_quad+1,1)];
                    f_ode = [v;invM*f_v;1];
                    % updated with clock state
                    F = [f_ode, f_aux_n1];
                    S = [1; -1];
                    n_aux = 1;
                    c = f_c;
                end

                %% Settings updates
                settings.time_freezing_model_exists = 1;
                %% Model updates
                model.dims.n_quad = n_quad;
                model.dims.n_q = n_q;
                model.dims.n_aux = n_aux;
                model.q = q;
                model.v = v;
                model.x = x;
                model.x0 = x0;
                model.M = M;
                model.dims.n_contacts = n_contacts;
                model.mu = mu;
                model.J_normal = J_normal;
                model.F = F;
                model.c = c;
                model.S = S;
                model.dims.n_dim_contact = 2;
            end
        end
    end  % methods
end % NosnocModel
