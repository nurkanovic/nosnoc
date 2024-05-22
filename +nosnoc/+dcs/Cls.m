classdef Cls < nosnoc.dcs.Base
    properties
        model

        % CLS stage vars
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
        z_v

        % CLS impulse vars
        Lambda_normal
        Y_gap
        L_vn
        P_vn
        N_vn
        Lambda_tangent
        Gamma_d
        Beta_d
        Delta_d
        Gamma
        Beta
        P_vt
        N_vt
        Alpha_vt

        g_lift

        z_alg
        z_impulse

        f_x

        M_fun
        invM_fun
        f_c_fun
        g_impulse_fun
        J_normal_fun
        J_tangent_fun
        D_tangent_fun

        dims
    end

    methods
        function obj = Cls(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            obj.lambda_normal = define_casadi_symbolic(opts.casadi_symbolic_mode,'lambda_normal',dims.n_c);
            obj.y_gap = define_casadi_symbolic(opts.casadi_symbolic_mode,'y_gap',dims.n_c);
            % Variables for impulse equations
            obj.Lambda_normal = define_casadi_symbolic(opts.casadi_symbolic_mode,'Lambda_normal',dims.n_c);
            obj.Y_gap = define_casadi_symbolic(opts.casadi_symbolic_mode,'Y_gap',dims.n_c);
            %obj.L_vn = define_casadi_symbolic(opts.casadi_symbolic_mode,'L_vn',dims.n_c); % lifting variable for state jump law
            obj.P_vn = define_casadi_symbolic(opts.casadi_symbolic_mode,'P_vn',dims.n_c); % lifting variable for state jump law
            obj.N_vn = define_casadi_symbolic(opts.casadi_symbolic_mode,'N_vn',dims.n_c); % lifting variable for state jump law
            if model.friction_exists
                % tangetial contact force (friction force)
                obj.lambda_tangent = define_casadi_symbolic(opts.casadi_symbolic_mode,'lambda_tangent',dims.n_tangents);
                % Impulse variables
                obj.Lambda_tangent = define_casadi_symbolic(opts.casadi_symbolic_mode,'Lambda_tangent',dims.n_tangents);
                if isequal(opts.friction_model,'Polyhedral')
                    obj.gamma_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'gamma_d',dims.n_c);
                    obj.beta_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'beta_d',dims.n_c); % lift friction cone bound
                    obj.delta_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'delta_d',dims.n_tangents); % lift lagrangian
                                                                                                               % Impulse variables
                    obj.Gamma_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'Gamma_d',dims.n_c);
                    obj.Beta_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'Beta_d',dims.n_c); % lift friction cone bound
                    obj.Delta_d = define_casadi_symbolic(opts.casadi_symbolic_mode,'Delta_d',dims.n_tangents); % lift lagrangian
                end
                if isequal(opts.friction_model,'Conic')
                    obj.gamma = define_casadi_symbolic(opts.casadi_symbolic_mode,'gamma',dims.n_c);
                    obj.beta = define_casadi_symbolic(opts.casadi_symbolic_mode,'beta',dims.n_c);
                    % Impulse variables
                    obj.Gamma = define_casadi_symbolic(opts.casadi_symbolic_mode,'Gamma',dims.n_c);
                    obj.Beta = define_casadi_symbolic(opts.casadi_symbolic_mode,'Beta',dims.n_c);
                    switch opts.conic_model_switch_handling
                      case 'Plain'
                        % no extra constraints
                      case 'Abs'
                        obj.p_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'p_vt',dims.n_tangents); % positive parts of tagnetial velocity (for switch detection)
                        obj.n_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'n_vt',dims.n_tangents); % negative parts of tagnetial velocity (for switch detection)
                                                                                                             % Impulse
                        obj.P_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'P_vt',dims.n_tangents);
                        obj.N_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'N_vt',dims.n_tangents);
                      case 'Lp'
                        obj.p_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'p_vt',dims.n_tangents); % positive parts of tagnetial velocity (for switch detection)
                        obj.n_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'n_vt',dims.n_tangents); % negative parts of tagnetial velocity (for switch detection)
                        obj.alpha_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'alpha_vt ',dims.n_tangents); % step function of tangential velocities
                                                                                                                      % impulse
                        obj.P_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'P_vt',dims.n_tangents);
                        obj.N_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'N_vt',dims.n_tangents);
                        obj.Alpha_vt = define_casadi_symbolic(opts.casadi_symbolic_mode,'Alpha_vt ',dims.n_tangents);
                    end
                end
            end
            if opts.lift_velocity_state
                obj.z_v = define_casadi_symbolic(casadi_symbolic_mode,['z_v'],dims.n_q);
            end
            
            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            if ~opts.lift_velocity_state
                if model.friction_exists
                    switch opts.friction_model
                      case 'Conic'
                        F_v = inv(model.M)*(model.f_v+model.J_normal*obj.lambda_normal+model.J_tangent*obj.lambda_tangent);
                      case 'Polyhedral'
                        F_v = inv(model.M)*(model.f_v+model.J_normal*obj.lambda_normal+model.D_tangent*obj.lambda_tangent);
                    end
                else
                    F_v = inv(model.M)*(model.f_v+model.J_normal*obj.lambda_normal);
                end
                obj.f_x = [model.v;F_v];
            else
                obj.f_x = [v;obj.z_v];
                if model.friction_exists
                    switch opts.friction_model
                      case 'Conic'
                        g_lift_v = model.M*obj.z_v -(model.f_v +model.J_normal*obj.lambda_normal + model.J_tangent*obj.lambda_tangent);
                      case 'Polyhedral'
                        g_lift_v =  model.M*obj.z_v -(model.f_v + model.J_normal*obj.lambda_normal + model.D_tangent*obj.lambda_tangent);
                    end
                else
                    g_lift_v =  model.M*obj.z_v -(model.f_v+model.J_normal*obj.lambda_n);
                end
                obj.g_lift = [obj.g_lift;g_lift_v];
            end
            
            lambda00_expr =[];
            
            % dummy variables for impact quations:
            v_post_impact = define_casadi_symbolic(opts.casadi_symbolic_mode,'v_post_impact',dims.n_q);
            v_pre_impact = define_casadi_symbolic(opts.casadi_symbolic_mode,'v_pre_impact',dims.n_q);
            g_alg_cls = [obj.y_gap - model.f_c];
            g_impulse = [model.M*(v_post_impact-v_pre_impact)-model.J_normal*obj.Lambda_normal]; % TODO: can this be relaxed? velocity junction
            g_impulse = [g_impulse; obj.Y_gap-model.f_c];
            % add state jump for every contact
            for ii = 1:dims.n_c
                %g_impulse = [g_impulse; obj.L_vn(ii) - model.J_normal(:,ii)'*(v_post_impact+model.e(ii)*v_pre_impact)];
                g_impulse = [g_impulse; obj.P_vn(ii) - obj.N_vn(ii) - model.J_normal(:,ii)'*(v_post_impact+model.e(ii)*v_pre_impact)];
            end
            if model.friction_exists
                switch opts.friction_model
                  case 'Conic'
                    g_impulse(1:dims.n_q) =  model.M*(v_post_impact-v_pre_impact)-model.J_normal*obj.Lambda_normal-model.J_tangent*obj.Lambda_tangent;
                    % algebraic and friction equations
                    for ii = 1:dims.n_c
                        ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                        g_alg_cls  = [g_alg_cls;-model.J_tangent(:,ind_temp)'*model.v - 2*obj.gamma(ii)*obj.lambda_tangent(ind_temp);...
                            obj.beta(ii)-((model.mu(ii)*obj.lambda_normal(ii))^2-norm(obj.lambda_tangent(ind_temp))^2)];
                        g_impulse = [g_impulse;
                            -model.J_tangent(:,ind_temp)'*v_post_impact - 2*obj.Gamma(ii)*obj.Lambda_tangent(ind_temp);...
                            obj.Beta - ((model.mu(ii)*obj.Lambda_normal(ii))^2- norm(obj.Lambda_tangent(ind_temp))^2)];

                        if ~isequal(opts.conic_model_switch_handling,'Plain')
                            % equality constraints for pos and neg parts of the tangetial velocity
                            g_alg_cls  = [g_alg_cls;model.J_tangent(:,ind_temp)'*model.v-(obj.p_vt(ind_temp)-obj.n_vt(ind_temp))];
                            g_impulse = [g_impulse;model.J_tangent(:,ind_temp)'*v_post_impact - (obj.P_vt(ind_temp)-obj.N_vt(ind_temp))];
                        end
                    end
                  case 'Polyhedral'
                    g_impulse(1:dims.n_q) = model.M*(v_post_impact-v_pre_impact)-model.J_normal*obj.Lambda_normal-model.D_tangent*obj.Lambda_tangent;
                    % impulse lifting equations
                    for ii = 1:dims.n_c
                        ind_temp = dims.n_t*ii-(dims.n_t-1):dims.n_t*ii;
                        g_alg_cls  = [g_alg_cls;obj.beta_d(ii)-(model.mu(ii)*obj.lambda_normal(ii) - sum(obj.lambda_tangent(ind_temp)));...
                            obj.delta_d(ind_temp) - (model.D_tangent(:,ind_temp)'*model.v + obj.gamma_d(ii))];
                        g_impulse = [g_impulse;obj.Beta_d - (model.mu(ii)*obj.Lambda_normal(ii)-sum(obj.Lambda_tangent(ind_temp)));...
                            obj.Delta_d(ind_temp)- (model.D_tangent(:,ind_temp)'*v_post_impact + obj.Gamma_d(ii))];
                    end
                end
            end
            g_alg = g_alg_cls;

            z_alg = [obj.lambda_normal; obj.y_gap; obj.lambda_tangent; obj.gamma_d; obj.beta_d; obj.delta_d; obj.beta_conic; obj.gamma_conic; obj.p_vt; obj.n_vt; obj.alpha_vt];
            %z_impulse = [obj.Lambda_normal; obj.Y_gap; obj.L_vn; obj.Lambda_tangent; obj.Gamma_d; obj.Beta_d; obj.Delta_d; obj.Gamma; obj.Beta; obj.P_vt; obj.N_vt; obj.Alpha_vt];
            z_impulse = [obj.Lambda_normal; obj.Y_gap; obj.P_vn; obj.N_vn; obj.Lambda_tangent; obj.Gamma_d; obj.Beta_d; obj.Delta_d; obj.Gamma; obj.Beta; obj.P_vt; obj.N_vt; obj.Alpha_vt];
            z_alg_f_x = [obj.lambda_normal; obj.lambda_tangent; obj.z_v];

            obj.f_x_fun = Function('f_x', {model.x, model.z, z_alg_f_x, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, z_alg, model.v_global, model.p}, {g_alg});
            obj.g_impulse_fun = Function('g_impulse', {model.q, v_post_impact, v_pre_impact, z_impulse, model.v_global, model.p}, {g_impulse}); % TODO (@anton) user algebraics pre and post impact?
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.g_comp_path_fun  = Function('g_comp_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_comp_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});

            obj.M_fun = Function('M_fun', {model.x}, {model.M});
            obj.invM_fun = Function('invM_fun', {model.x}, {model.invM});
            obj.f_c_fun = Function('f_c_fun', {model.x}, {model.f_c});
            obj.J_normal_fun = Function('J_normal_fun', {model.x}, {model.J_normal});
            if model.friction_exists
                if isequal(opts.friction_model,'Conic')
                    obj.J_tangent_fun = Function('J_tangent_fun', {model.x}, {model.J_tangent});
                else
                    obj.D_tangent_fun = Function('D_tangent_fun', {model.x}, {model.D_tangent});
                end
            end
        end
    end
end
