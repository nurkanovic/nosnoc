classdef Heaviside < nosnoc.dcs.Base
    properties
        alpha % CasADi symbolic variable for selection of the Heaviside step function
        alpha_sys % cell containing the alpha variables of every subsystem, wheras alpha stores the concatenation of all these vectors;
        lambda_n % CasADi symbolic variable 
        lambda_n_sys % cell
        lambda_p % CasADi symbolic variable 
        lambda_p_sys % cell
        % These are relevant only for lifting. For now wait to implement until re-implementing time-freezing
        % beta 
        % gamma
        % theta
        % theta_sys
        theta_expr_sys 

        z_all

        f_x
        
        dims

        g_lp_stationarity_fun
        lambda00_fun
    end

    methods
        function obj = Heaviside(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            if class(obj.model) == "nosnoc.model.Cls"
                obj.time_freezing(opts);
            end
            dims = obj.dims;
            model = obj.model;

            if ~ismember(class(obj.model), {'nosnoc.model.Pss', 'nosnoc.model.Heaviside'})
                error('A heaviside DCS can only be created from a PSS or Heaviside step model.')
            end

            switch class(obj.model)
              case "nosnoc.model.Pss"
                % Check pss popluated S and c;
                if isempty(model.S)
                    error("The PSS model must containt the lp_stationarity matrix S and lp_stationarity function c.");
                end

                dims.n_alpha = sum(obj.dims.n_c_sys); % number of modes
              case "nosnoc.model.Heaviside"
                obj.alpha = model.alpha;
            end
            dims.n_lambda = dims.n_alpha;

            idx = 1;
            for ii = 1:dims.n_sys
                sys_idx{ii} = idx:(idx + obj.dims.n_c_sys(ii)-1);
                idx = idx + obj.dims.n_c_sys(ii);
                ii_str = num2str(ii);
                % define alpha (selection of a set valued step function)
                switch class(obj.model)
                  case "nosnoc.model.Pss"
                    obj.alpha_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                    obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                  case "nosnoc.model.Heaviside"
                    obj.alpha_sys{ii} = obj.alpha(sys_idx{ii});
                end
                % define lambda_n_i (Lagrange multipler of alpha >= 0;)
                obj.lambda_n_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['lambda_n_' ii_str],obj.dims.n_c_sys(ii));
                obj.lambda_n = [obj.lambda_n;obj.lambda_n_sys{ii}];
                % define lambda_p_i (Lagrange multipler of alpha <= 1;)
                obj.lambda_p_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['lambda_p_' ii_str],obj.dims.n_c_sys(ii));
                obj.lambda_p = [obj.lambda_p;obj.lambda_p_sys{ii}];
            end

            % symbolic variables z = [alpha;obj.lambda_n;obj.lambda_p];
            obj.z_all = [obj.alpha;obj.lambda_n;obj.lambda_p;obj.model.z];
            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            switch class(obj.model)
              case "nosnoc.model.Pss"
                % generate theta_expressions
                for ii = 1:dims.n_sys
                    theta_i = [];
                    S_temp = model.S{ii};
                    if opts.pss_lift_step_functions
                        % TODO implement automatic lifting
                    else
                        for jj = 1:size(S_temp,1)
                            theta_jk = 1;
                            for kk = 1:size(S_temp,2)
                                % create multiaffine term
                                if S_temp(jj,kk) ~=0
                                    theta_jk = theta_jk * (0.5*(1-S_temp(jj,kk))+S_temp(jj,kk)*obj.alpha_sys{ii}(kk));
                                end
                            end
                            theta_i = [theta_i;theta_jk];
                        end                    
                    end
                    obj.theta_expr_sys{ii} = theta_i;
                end

                obj.f_x = zeros(dims.n_x,1);
                for ii = 1:dims.n_sys
                    obj.f_x = obj.f_x + model.F{ii}*obj.theta_expr_sys{ii};
                end
              case "nosnoc.model.Heaviside"
                obj.f_x = model.f_x;
            end

            g_lp_stationarity = []; % collects lp_stationarity function algebraic equations, 0 = c(x)-lambda_p+lambda_n
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                % c_i(x) - (lambda_p_i-lambda_n_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i'*alpha_i  = 0; for all i = 1,..., n_sys
                % lambda_p_i'*(e-alpha_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i >= 0;    for all i = 1,..., n_sys
                % lambda_p_i >= 0;    for all i = 1,..., n_sys
                % alpha_i >= 0;     for all i = 1,..., n_sys
                g_lp_stationarity = [g_lp_stationarity;model.c{ii}-obj.lambda_p_sys{ii}+obj.lambda_n_sys{ii}];
                lambda00_expr = [lambda00_expr; -min(model.c{ii}, 0); max(model.c{ii},0)];
            end
            g_alg = [g_lp_stationarity];

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {g_alg});
            obj.g_lp_stationarity_fun = Function('g_lp_stationarity', {model.x, model.z, obj.lambda_n, obj.lambda_p, model.v_global, model.p}, {g_lp_stationarity});
            obj.lambda00_fun = Function('lambda00', {model.x, model.z, model.v_global, model.p_global}, {lambda00_expr});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.G_path_fun  = Function('G_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.G_path});
            obj.H_path_fun  = Function('H_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.H_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});
        end

        function time_freezing(obj, opts)
            model = obj.model;
            dims = obj.dims;
            pss_model = nosnoc.model.Pss();
            pss_model.dims = obj.dims;

            % transfer all common properties to pss_model
            common_props = properties('nosnoc.model.Base')
            for ii=1:length(common_props)
                prop = common_props{ii};
                pss_model.(prop) = model.(prop);
            end
            
            dims.n_contacts = length(model.f_c);
            % check dimensions of contacts
            if isempty(dims.n_dim_contact)
                warning('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
                dims.n_dim_contact = 2;
            end
            
            % qudrature state
            % TODO remember to add to derivative in the end
            dims.n_quad = 0;
            f_quad = [];
            if opts.time_freezing_quadrature_state
                % define quadrature state
                L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
                pss_model.lbx = [pss_model.lbx;-inf];
                pss_model.ubx = [pss_model.ubx;inf];
                pss_model.x = [pss_model.x;L];
                pss_model.x0 = [pss_model.x0;0];
                f_quad = pss_model.f_q;
                pss_model.f_q = 0;
                if ~isempty(pss_model.f_q_T)
                    pss_model.f_q_T  = pss_model.f_q_T + L;
                else
                    pss_model.f_q_T = L;
                end
                dims.n_quad = 1;
            end
          
            % uneven number of states = it is assumed that the clock state is defined.
            t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
            % update lower and upper bounds of lbx and ubx
            pss_model.lbx = [pss_model.lbx;-inf];
            pss_model.ubx = [pss_model.ubx;inf];
            pss_model.x = [pss_model.x;t];
            pss_model.x0 = [pss_model.x0;0];
            
            % normal and tangential velocities
            eps_t = 1e-7;
            v_normal = model.J_normal'*model.v;
            if obj.friction_exists
                if dims.n_dim_contact == 2
                    v_tangent = (model.J_tangent'*model.v)';
                else
                    v_tangent = model.J_tangent'*model.v;
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
            % TODO add to options
            a_n = 100;
            %% Time-freezing reformulation
            if model.e == 0
                % Basic problem_options
                %% switching function
                if opts.nonsmooth_switching_fun
                    pss_model.c = [max_smooth_fun(model.f_c,v_normal,0);v_tangent];
                else
                    if dims.n_dim_contact == 2
                        pss_model.c = [model.f_c;v_normal;v_tangent'];
                    else
                        pss_model.c = [model.f_c;v_normal;v_tangent_norms-eps_t];
                    end
                end
                %% unconstrained dynamics with clock state
                f_ode = [model.v;...
                    model.invM*model.f_v;
                    f_quad;
                    1];

                %% Auxiliary dynamics
                % where to use invM, in every auxiliary dynamics or only at the end
                if inv_M_once
                    inv_M_aux = eye(dims.n_q);
                    inv_M_ext = blkdiag(zeros(dims.n_q),model.invM,0);
                else
                    inv_M_aux = model.invM;
                    inv_M_ext = eye(dims.n_x+1);
                end
                f_aux_pos = []; % matrix with all auxiliary tangent dynamics
                f_aux_neg = [];
                % time freezing dynamics
                if opts.stabilizing_q_dynamics
                    f_q_dynamics = -opts.kappa_stabilizing_q_dynamics*obj.J_normal*diag(model.f_c);
                else
                    f_q_dynamics = zeros(dims.n_q,dims.n_contacts);
                end
                f_aux_normal = [f_q_dynamics;inv_M_aux*obj.J_normal*obj.a_n;zeros(1,dims.n_contacts)];

                for ii = 1:dims.n_contacts
                    if model.friction_exists && model.mu_f(ii)>0
                        % auxiliary tangent;
                        if dims.n_dim_contact == 2
                            v_tangent_ii = model.J_tangent(:,ii)'*model.v;
                            f_aux_pos_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(model.J_normal(:,ii)-model.J_tangent(:,ii)*(model.mu_f(ii)))*a_n;0]; % for v>0
                            f_aux_neg_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(model.J_normal(:,ii)+model.J_tangent(:,ii)*(model.mu_f(ii)))*a_n;0]; % for v<0
                        else
                            v_tangent_ii = v_tangent(:,ii);
                            f_aux_pos_ii = [f_q_dynamics(:,ii);inv_M_aux*(model.J_normal(:,ii)*a_n-model.J_tangent(:,ii*2-1:ii*2)*model.mu_f(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                            f_aux_neg_ii = [f_q_dynamics(:,ii);inv_M_aux*(model.J_normal(:,ii)*a_n+model.J_tangent(:,ii*2-1:ii*2)*model.mu_f(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                        end
                        f_aux_pos = [f_aux_pos,f_aux_pos_ii];
                        f_aux_neg = [f_aux_neg,f_aux_neg_ii];
                    end
                end
                if model.friction_exists
                    f_aux = [f_aux_pos,f_aux_neg];
                else
                    f_aux = f_aux_normal;
                end
                pss_model.F = [f_ode (inv_M_ext*f_aux)];
                pss.S = ones(size(obj.F,2),length(obj.c)); % dummy value to pass error checks
                                                           % number of auxiliary dynamicsm modes
                if model.friction_exists
                    dims.n_aux = 2*dims.n_contacts;
                else
                    dims.n_aux = dims.n_contacts;
                end
            else
                % elastic
                % TODO: add to options
                if true
                    k_aux = 10;
                    if opts.print_level > 1
                        fprintf('nosnoc: Setting default value for k_aux = 10.\n')
                    end
                end
                temp1 = 2*abs(log(model.e));
                temp2 = k_aux/(pi^2+log(model.e)^2);
                c_aux = temp1/sqrt(temp2);
                K = [0 1;-k_aux -c_aux];
                N  = [obj.J_normal zeros(dims.n_q,1);...
                    zeros(dims.n_q,1) obj.invM*obj.J_normal];
                f_aux_n1 = N*K*N'*[model.q;model.v];
                f_aux_n1 = [f_aux_n1;zeros(dims.n_quad+1,1)];
                f_ode = [model.v;model.invM*model.f_v;f_quad;1];
                % updated with clock state
                pss_model.F = [f_ode, f_aux_n1];
                pss_model.S = [1; -1];
                pss_model.c = obj.f_c;
                dims.n_aux = 1;
            end

            obj.dims = dims;
        end
    end

    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            propgrp = getPropertyGroups@nosnoc.dcs.Base(obj);
            group_title = 'Variables';
            var_list = struct;
            var_list.lambda_p = obj.lambda_p;
            var_list.lambda_n = obj.lambda_n;
            var_list.alpha = obj.alpha;
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end

end
