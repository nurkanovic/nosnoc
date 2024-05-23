classdef Cls < vdx.problems.Mpcc
    properties
        model
        dcs
        opts

        populated
    end

    properties(Access=private)
        z_alg
        z_impulse
        z_alg_f_x
    end

    methods
        function obj = Cls(dcs, opts)
            obj = obj@vdx.problems.Mpcc();
            obj.model = dcs.model;
            obj.dcs = dcs;
            obj.opts = opts;
        end

        function create_variables(obj)
            dims = obj.dcs.dims;
            dcs = obj.dcs;
            model = obj.model;
            opts = obj.opts;

            % Parameters
            obj.p.rho_h_p = {{'rho_h_p',1}, 1};
            obj.p.rho_terminal_p = {{'rho_terminal_p',1}, 1};
            obj.p.T = {{'T',1}, opts.T};
            obj.p.p_global = {model.p_global, model.p_global_val};


            % 0d vars: Variables which only exist once globally. 
            obj.w.v_global = {{'v_global',dims.n_v_global}, model.lbv_global, model.ubv_global, model.v0_global};
            if opts.use_speed_of_time_variables && ~opts.local_speed_of_time_variable
                obj.w.sot = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
            end
            if opts.time_optimal_problem
                obj.w.T_final = {{'T_final', 1}, opts.T_final_min, opts.T_final_max, opts.T};
                obj.f = obj.f + obj.w.T_final();
            end

            % 1d vars: Variables that are defined per control stage
            obj.w.u(1:opts.N_stages) = {{'u', dims.n_u}, model.lbu, model.ubu, model.u0};
            obj.p.p_time_var(1:opts.N_stages) = {{'p_time_var', dims.n_p_time_var}, model.p_time_var_val};
            if opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                obj.w.sot(1:opts.N_stages) = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
            end

            % TODO(@anton) This _severely_ hurts performance over the vectorized assignment by doing N_stages vertcats of
            %              casadi symbolics vs just a vectorized assignment which does one. As such there needs to be backend
            %              work done for vdx to cache vertcats of SX somehow. Current theory is one can simply keep a queue of
            %              symbolics to be added in a cell array until a read is done, at which point we call a single vertcat
            %              on the whole queue which is _significantly_ faster.
            % 2d vars: Variables that are defined for each finite element.
            for ii=1:opts.N_stages
                % other derived values
                h0 = opts.h_k(ii);
                if obj.opts.use_fesd
                    ubh = (1 + opts.gamma_h) * h0;
                    lbh = (1 - opts.gamma_h) * h0;
                    if opts.time_rescaling && ~opts.use_speed_of_time_variables
                        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                        ubh = (1+opts.gamma_h)*h0*opts.s_sot_max;
                        lbh = (1-opts.gamma_h)*h0/opts.s_sot_min;
                    end
                    obj.w.h(ii,1:opts.N_finite_elements(ii)) = {{'h', 1}, lbh, ubh, h0};
                end
                if obj.opts.step_equilibration == StepEquilibrationMode.mlcp
                    % TODO fix this :)
                    obj.w.B_max(ii,2:opts.N_finite_elements(ii)) = {{'B_max', dims.n_lambda},-inf,inf};
                    obj.w.pi_c(ii,2:opts.N_finite_elements(ii)) = {{'pi_c', dims.n_c},-inf,inf};
                    obj.w.pi_lambda(ii,2:opts.N_finite_elements(ii)) = {{'pi_lambda', dims.n_lambda},-inf,inf};
                    obj.w.lambda_c(ii,2:opts.N_finite_elements(ii)) = {{'lambda_c', dims.n_c},0,inf};
                    obj.w.lambda_lambda(ii,2:opts.N_finite_elements(ii)) = {{'lambda_lambda', dims.n_lambda},0,inf};
                    obj.w.eta(ii,2:opts.N_finite_elements(ii)) = {{'eta', dims.n_lambda},0,inf};
                    obj.w.nu(ii,2:opts.N_finite_elements(ii)) = {{'nu', 1},0,inf};
                end
            end

            if opts.no_initial_impacts
                start_fe = 2;
            else
                start_fe = 1;
            end
            % 2d Impulse vars
            obj.w.Lambda_normal(ii,start_fe:opts.N_finite_elements(ii)) = {{'Lambda_normal', dims.n_c}, 0, inf, 1};
            %obj.w.L_vn(ii,start_fe:opts.N_finite_elements(ii)) = {{'L_vn', dims.n_c}, -inf, inf, 1};
            obj.w.P_vn(ii,start_fe:opts.N_finite_elements(ii)) = {{'P_vn', dims.n_c}, 0, inf, 1};
            obj.w.N_vn(ii,start_fe:opts.N_finite_elements(ii)) = {{'N_vn', dims.n_c}, 0, inf, 1};
            obj.w.Y_gap(ii,start_fe:opts.N_finite_elements(ii)) = {{'Y_gap', dims.n_c}, 0, inf, 1};
            if model.friction_exists
                switch opts.friction_model
                  case 'Polyhedral'
                    obj.w.Lambda_tangent(ii,start_fe:opts.N_finite_elements(ii)) = {{'Lambda_tangent', dims.n_tangents}, 0, inf, 1};
                    obj.w.Gamma_d(ii,start_fe:opts.N_finite_elements(ii)) = {{'Gamma_d', dims.n_c}, 0, inf, 1};
                    obj.w.Beta_d(ii,start_fe:opts.N_finite_elements(ii)) = {{'Beta_d', dims.n_c}, 0, inf, 1};
                    obj.w.Delta_d(ii,start_fe:opts.N_finite_elements(ii)) = {{'Delta_d', dims.n_tangents}, 0, inf, 1};
                  case 'Conic'
                    obj.w.Lambda_tangent(ii,start_fe:opts.N_finite_elements(ii)) = {{'Lambda_tangent', dims.n_tangents}, -inf, inf, 0};
                    obj.w.Gamma(ii,start_fe:opts.N_finite_elements(ii)) = {{'Gamma', dims.n_tangents}, 0, inf, 1};
                    obj.w.Beta(ii,start_fe:opts.N_finite_elements(ii)) = {{'Beta', dims.n_tangents}, 0, inf, 1};
                    switch opts.conic_model_switch_handling
                      case 'Plain'
                        % no extra vars
                      case 'Abs'
                        obj.w.P_vt(ii,start_fe:opts.N_finite_elements(ii)) = {{'P_vt', dims.n_tangents}, 0, inf, 1};
                        obj.w.N_vt(ii,start_fe:opts.N_finite_elements(ii)) = {{'N_vt', dims.n_tangents}, 0, inf, 1};
                      case 'Lp'
                        obj.w.P_vt(ii,start_fe:opts.N_finite_elements(ii)) = {{'P_vt', dims.n_tangents}, 0, inf, 1};
                        obj.w.N_vt(ii,start_fe:opts.N_finite_elements(ii)) = {{'N_vt', dims.n_tangents}, 0, inf, 1};
                        obj.w.Alpha_vt(ii,start_fe:opts.N_finite_elements(ii)) = {{'Alpha_vt', dims.n_tangents}, 0, inf, 1};
                    end
                end
            end

            % For c_n ~= 1 case
            rbp = ~opts.right_boundary_point_explicit;
            
            % 3d vars: Variables defined on each rk stage
            %          some of which are also defined at the initial point:
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}, model.x0, model.x0, model.x0};
            obj.w.z(0,0,opts.n_s) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
            if model.friction_exists
                switch opts.friction_model
                  case 'Polyhedral'
                    obj.w.gamma_d(0,0,opts.n_s) = {{'gamma_d', dims.n_c}, 0, inf, 0};
                    obj.w.delta_d(0,0,opts.n_s) = {{'delta_d', dims.n_tangents}, 0, inf, 0};
                  case 'Conic'
                    obj.w.gamma(0,0,opts.n_s) = {{'gamma', dims.n_tangents}, 0, inf, 0};
                    switch opts.conic_model_switch_handling
                      case 'Plain'
                        % nothing
                      case {'Abs', 'Lp'}
                        obj.w.p_vt(0,0,opts.n_s) = {{'p_vt', dims.n_tangents}, 0, inf, 0};
                        obj.w.n_vt(0,0,opts.n_s) = {{'n_vt', dims.n_tangents}, 0, inf, 0};
                    end
                end
            end            
            for ii=1:opts.N_stages
                if (opts.rk_representation == RKRepresentation.integral ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    obj.w.x(ii,1:opts.N_finite_elements(ii),0:(opts.n_s+rbp)) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                else
                    obj.w.x(ii,1:opts.N_finite_elements(ii),0) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                    obj.w.x(ii,1:opts.N_finite_elements(ii),opts.n_s) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                end
                if (opts.rk_representation == RKRepresentation.differential ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    obj.w.v(ii,1:opts.N_finite_elements(ii),1:opts.n_s) = {{'v', dims.n_x}};
                end

                if opts.lift_velocity_state
                    obj.w.z_v(ii,1:opts.N_finite_elements(ii),1:opts.n_s) = {{'z_v', dims.n_x}};
                end
                
                obj.w.z(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
                obj.w.lambda_normal(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'lambda_normal', dims.n_c}, 0, 50, 0};
                obj.w.y_gap(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'y_gap', dims.n_c}, 0, inf, 0};
                if rbp
                    obj.w.y_gap(ii,opts.N_finite_elements(ii),(opts.n_s+rbp)) = {{'y_gap', dims.n_c}, 0, inf, 0};
                end
                if model.friction_exists
                    switch opts.friction_model
                      case 'Polyhedral'
                        obj.w.lambda_tangent(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'lambda_tangent', dims.n_tangents}, 0, 50, 0};
                        obj.w.gamma_d(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'gamma_d', dims.n_c}, 0, inf, 0};
                        obj.w.beta_d(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'beta_d', dims.n_c}, 0, inf, 1};
                        obj.w.delta_d(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'delta_d', dims.n_tangents}, 0, inf, 1};
                      case 'Conic'
                        obj.w.lambda_tangent(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'lambda_tangent', dims.n_tangents}, -20, 20, 0};
                        obj.w.gamma(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'gamma', dims.n_tangents}, 0, inf, 1};
                        obj.w.beta(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'beta', dims.n_tangents}, 0, inf, 1};
                        switch opts.conic_model_switch_handling
                          case 'Plain'
                            % no extra vars
                          case 'Abs'
                            obj.w.p_vt(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'p_vt', dims.n_tangents}, 0, inf, 1};
                            obj.w.n_vt(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'n_vt', dims.n_tangents}, 0, inf, 1};
                          case 'Lp'
                            obj.w.p_vt(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'p_vt', dims.n_tangents}, 0, inf, 1};
                            obj.w.n_vt(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'n_vt', dims.n_tangents}, 0, inf, 1};
                            obj.w.alpha_vt(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'alpha_vt', dims.n_tangents}, 0, inf, 1};
                        end
                    end
                end
                
                % Handle x_box settings
                if ~opts.x_box_at_stg && opts.rk_representation ~= RKRepresentation.differential
                    obj.w.x(1:opts.N_stages,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(1:opts.N_stages,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).ub = inf*ones(dims.n_x, 1);
                end

                if ~opts.x_box_at_fe
                    obj.w.x(1:opts.N_stages,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(1:opts.N_stages,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).ub = inf*ones(dims.n_x, 1);
                end
            end

            % build z_alg, z_impulse, and z_alg_f_x variable groups
            z_alg = {obj.w.lambda_normal, obj.w.y_gap};
            z_impulse = {obj.w.Lambda_normal, obj.w.Y_gap, obj.w.P_vn, obj.w.N_vn};
            z_alg_f_x = {obj.w.lambda_normal};
            if model.friction_exists
                z_alg = [z_alg, {obj.w.lambda_tangent}];
                z_impulse = [z_impulse, {obj.w.Lambda_tangent}];
                z_alg_f_x = [z_alg_f_x, {obj.w.lambda_tangent}];
                switch opts.friction_model
                  case 'Polyhedral'
                    z_alg = [z_alg,{obj.w.gamma_d,obj.w.beta_d,obj.w.delta_d}];
                    z_impulse = [z_impulse,{obj.w.Gamma_d,obj.w.Beta_d,obj.w.Delta_d}];
                  case 'Conic'
                    z_alg = [z_alg,{obj.w.gamma,obj.w.beta}];
                    z_impulse = [z_impulse,{obj.w.Gamma,obj.w.Beta}];
                    switch opts.conic_model_switch_handling
                      case 'Plain'
                        % no extra vars
                      case 'Abs'
                        z_alg = [z_alg,{obj.w.p_vt,obj.w.n_vt}];
                        z_impulse = [z_impulse,{obj.w.P_vt,obj.w.N_vt}];
                      case 'Lp'
                        z_alg = [z_alg,{obj.w.p_vt,obj.w.n_vt,obj.w.alpha_vt}];
                        z_impulse = [z_impulse,{obj.w.P_vt,obj.w.N_vt,obj.w.Alpha_vt}];
                    end
                end
            end
            if opts.lift_velocity_state
                z_alg_f_x = [z_alg_f_x, {obj.w.z_v}]
            end
            obj.z_alg = vdx.VariableGroup(z_alg);
            obj.z_impulse = vdx.VariableGroup(z_impulse);
            obj.z_alg_f_x = vdx.VariableGroup(z_alg_f_x);
        end

        function generate_direct_transcription_constraints(obj)
            import casadi.*
            model = obj.model;
            opts = obj.opts;
            dcs = obj.dcs;
            dims = obj.dcs.dims;

            rbp = ~opts.right_boundary_point_explicit;

            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            x_0 = obj.w.x(0,0,opts.n_s);
            z_0 = obj.w.z(0,0,opts.n_s);
            
            obj.g.z(0,0,opts.n_s) = {dcs.g_z_fun(x_0, z_0, obj.w.u(1), v_global, [p_global;obj.p.p_time_var(1)])};
            
            x_prev = obj.w.x(0,0,opts.n_s);
            for ii=1:opts.N_stages
                if obj.opts.use_fesd
                    t_stage = obj.p.T()/opts.N_stages;
                    h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(ii));
                else
                    h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(ii));
                end
                
                ui = obj.w.u(ii);
                p_stage = obj.p.p_time_var(ii);
                p = [p_global;p_stage];
                if obj.opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                    s_sot = obj.w.sot(ii);
                elseif obj.opts.use_speed_of_time_variables
                    s_sot = obj.w.sot();
                else
                    s_sot = 1;
                end

                sum_h = 0;
                for jj=1:opts.N_finite_elements(ii)
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    elseif opts.time_optimal_problem && ~opts.use_speed_of_time_variables
                        h = obj.w.T_final()/(opts.N_stages*opts.N_finite_elements(ii));
                    else
                        h = h0;
                    end

                    % left boundary point handling
                    q_prev = x_prev(1:dims.n_q);
                    v_prev = x_prev(dims.n_q + 1:end);
                    x_lbp = obj.w.x(ii,jj,0);
                    q_lbp = x_lbp(1:dims.n_q);
                    v_lbp = x_lbp(dims.n_q + 1:end);
                    obj.g.q_continuity(ii,jj) = {q_prev - q_lbp};
                    z_impulse_ij = obj.z_impulse(ii,jj);

                    if (jj ~= 1 || ~opts.no_initial_impacts)
                        obj.g.impulse(ii,jj) = {dcs.g_impulse_fun(q_lbp,v_lbp,v_prev,z_impulse_ij, v_global, p)};
                    else
                        obj.g.v_continuity(ii,jj) = {v_lbp-v_prev};
                    end
                    % additional y_eps constraint
                    if opts.eps_cls > 0
                        x_eps = vertcat(q_lbp + h * opts.eps_cls * v_lbp, v_lbp);
                        obj.g.f_c_eps(ii,jj) = {dcs.f_c_fun(x_eps), 0, inf};
                    end
                    
                    switch opts.rk_representation
                      case RKRepresentation.integral
                        % In integral representation stage variables are states.
                        x_ij_end = x_lbp;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            z_alg_ijk = obj.z_alg(ii,jj,kk);
                            z_alg_f_x_ijk = obj.z_alg_f_x(ii,jj,kk);

                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, z_alg_f_x_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, ui, v_global, p);
                            xk = opts.C_rk(1, kk+1) * x_lbp;
                            for rr=1:opts.n_s
                                x_ijr = obj.w.x(ii,jj,rr);
                                xk = xk + opts.C_rk(rr+1, kk+1) * x_ijr;
                            end
                            obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, z_alg_ijk, v_global, p)};

                            x_ij_end = x_ij_end + opts.D_rk(kk+1)*x_ijk;
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.g_comp_path, 1) > 0
                                g_comp_path = dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {g_comp_path(:,1)};
                                obj.H.path(ii,jj,kk) = {g_comp_path(:,2)};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                % TODO(@anton) make an option
                                %obj.f = obj.f + opts.B_rk(kk+1)*h*(qj + 1e-7*obj.w.lambda_normal(ii,jj,kk).^2);
                                obj.f = obj.f + opts.B_rk(kk+1)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            y_gap_ijk = obj.w.z(ii,jj,opts.n_s+1);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.y_gap_rbp(ii.jj) = {y_gap_ijk - dcs.f_c_fun(x_ijk)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                        end
                      case RKRepresentation.differential
                        % In differential representation stage variables are the state derivatives.
                        X_ijk = {};
                        for kk = 1:opts.n_s
                            x_temp = x_lbp;
                            for rr = 1:opts.n_s
                                x_temp = x_temp + h*opts.A_rk(kk,rr)*obj.w.v(ii,jj,rr);
                            end
                            X_ijk = [X_ijk {x_temp}];
                        end
                        X_ijk = [X_ijk, {obj.w.x(ii,jj,opts.n_s)}];
                        x_ij_end = x_lbp;
                        for kk=1:opts.n_s
                            x_ijk = X_ijk{kk};
                            v_ijk = obj.w.v(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            z_alg_ijk = obj.z_alg(ii,jj,kk);
                            z_alg_f_x_ijk = obj.z_alg_f_x(ii,jj,kk);
                            
                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, z_alg_f_x_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, ui, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {fj - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, z_alg_ijk, v_global, p)};
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.g_comp_path, 1) > 0
                                g_comp_path = dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {g_comp_path(:,1)};
                                obj.H.path(ii,jj,kk) = {g_comp_path(:,2)};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            y_gap_ijk = obj.w.y_gap(ii,jj,opts.n_s+1);
                            
                            
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.y_gap_rbp(ii.jj) = {y_gap_ijk - dcs.f_c_fun(x_ijk)};
                        else
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ij_end - obj.w.x(ii,jj,opts.n_s)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                        end
                      case RKRepresentation.differential_lift_x
                        % In differential representation with lifted state stage variables are the state derivatives and we
                        % lift the states at each stage point as well.
                        for kk = 1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            x_temp = x_lbp;
                            for rr = 1:opts.n_s
                                x_temp = x_temp + h*opts.A_rk(kk,rr)*obj.w.v(ii,jj,rr);
                            end
                            obj.g.lift_x(ii,jj,kk) = {x_ijk - x_temp};
                        end
                        x_ij_end = x_lbp;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            v_ijk = obj.w.v(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            z_alg_ijk = obj.z_alg(ii,jj,kk);
                            z_alg_f_x_ijk = obj.z_alg_f_x(ii,jj,kk);
                           
                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, z_alg_f_x_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, ui, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {fj - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, z_alg_ijk, v_global, p)};
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.g_comp_path, 1) > 0
                                g_comp_path = dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {g_comp_path(:,1)};
                                obj.H.path(ii,jj,kk) = {g_comp_path(:,2)};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            y_gap_ijk = obj.w.y_gap(ii,jj,opts.n_s+1);
                            
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.y_gap_rbp(ii.jj) = {y_gap_ijk - dcs.f_c_fun(x_ijk)};
                        else
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ij_end - obj.w.x(ii,jj,opts.n_s)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                        end
                    end
                    x_prev = obj.w.x(ii,jj,opts.n_s+rbp);
                end
                if ~opts.g_path_at_stg && ~opts.g_path_at_fe
                    x_i = obj.w.x(ii, opts.N_finite_elements(ii), opts.n_s);
                    z_i = obj.w.z(ii, opts.N_finite_elements(ii), opts.n_s);
                    obj.g.path(ii) = {dcs.g_path_fun(x_i, z_i, ui, v_global, p), model.lbg_path, model.ubg_path};
                end

                % Least Squares Costs
                % TODO we should convert the refs to params
                obj.f = obj.f + t_stage*dcs.f_lsq_x_fun(obj.w.x(ii,opts.N_finite_elements(ii),opts.n_s),...
                    model.x_ref_val(:,ii),...
                    p);
                obj.f = obj.f + t_stage*dcs.f_lsq_u_fun(obj.w.u(ii),...
                    model.u_ref_val(:,ii),...
                    p);
                if ~opts.cost_integration
                    obj.f = obj.f + dcs.f_q_fun(x_ijk, z_ijk, ui, v_global, p);
                end

                % Clock <Constraints
                % TODO(@anton) HERE BE DRAGONS. This is by far the worst part of current nosnoc as it requires the discrete problem
                %              to understand something about the time-freezing reformulation which is ugly.
                if obj.opts.use_fesd && opts.equidistant_control_grid
                    if opts.time_optimal_problem
                        if opts.use_speed_of_time_variables
                            obj.g.equidistant_control_grid(ii) = {[sum_h - opts.h;sot*sum_h - obj.w.T_final()/opts.N_stages]};
                        else
                            obj.g.equidistant_control_grid(ii) = {sum_h - obj.w.T_final()/opts.N_stages};
                        end
                    else
                        if opts.relax_terminal_numerical_time
                            % TODO(@armin) why is this allowed to be negative?
                            obj.w.s_numerical_time(ii) = {{'s_numerical', 1}, -2*opts.h, 2*opts.h, opts.h/2};
                            g_eq_grid = [sum_h - t_stage - obj.w.s_numerical_time(ii);
                                -(sum_h - t_stage) - obj.w.s_numerical_time(ii)];
                            obj.g.equidistant_control_grid(ii) = {g_eq_grid, -inf, 0};
                            obj.f = obj.f + opts.rho_terminal_numerical_time*obj.w.s_numerical_time(ii);
                        else
                            obj.g.equidistant_control_grid(ii) = {t_stage-sum_h};
                        end
                    end
                end
            end

            x_end = obj.w.x(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            z_end = obj.w.z(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            % Terminal cost
            obj.f = obj.f + dcs.f_q_T_fun(x_end, z_end, v_global, p_global);

            % Terminal_lsq_cost
            obj.f = obj.f + h0*opts.N_finite_elements(ii)*dcs.f_lsq_T_fun(x_end,...
                model.x_ref_end_val,...
                p_global);

            % Terminal constraint
            if opts.relax_terminal_constraint_homotopy
                error("Currently unsupported")
            end
            g_terminal = dcs.g_terminal_fun(x_end, z_end, v_global, p_global);
            switch opts.relax_terminal_constraint
              case ConstraintRelaxationMode.NONE % hard constraint
                if opts.relax_terminal_constraint_from_above
                    obj.g.terminal = {g_terminal, model.lbg_terminal, inf*ones(dims.n_g_terminal,1)};
                else
                    obj.g.terminal = {g_terminal, model.lbg_terminal, model.ubg_terminal};
                end
              case ConstraintRelaxationMode.ELL_1 % l_1
                obj.w.s_terminal_ell_1 = {{'s_terminal_ell_1', dims.n_g_terminal}, 0, inf, 10};

                g_terminal = [g_terminal-model.lbg_terminal-obj.w.s_terminal_ell_1();
                    -(g_terminal-model.ubg_terminal)-obj.w.s_terminal_ell_1()];
                obj.g.terminal = {g_terminal, -inf, 0};
                obj.f = obj.f + obj.p.rho_terminal_p()*sum(obj.w.s_terminal_ell_1());
              case ConstraintRelaxationMode.ELL_2 % l_2
                     % TODO(@anton): this is as it was implemented before. should handle lb != ub?
                obj.f = obj.f + obj.p.rho_terminal_p()*(g_terminal-model.lbg_terminal)'*(g_terminal-model.lbg_terminal);
              case ConstraintRelaxationMode.ELL_INF % l_inf
                obj.w.s_terminal_ell_inf = {{'s_terminal_ell_inf', 1}, 0, inf, 1e3};

                g_terminal = [g_terminal-model.lbg_terminal-obj.w.s_terminal_ell_inf();
                    -(g_terminal-model.ubg_terminal)-obj.w.s_terminal_ell_inf()];
                obj.g.terminal = {g_terminal, -inf, 0};
                obj.f = obj.f + obj.p.rho_terminal_p()*obj.w.s_terminal_ell_inf();
            end
        end

        % TODO(@anton) is there a way to make this less horribly ugly?
        function generate_complementarity_constraints(obj)
            import casadi.*
            opts = obj.opts;
            dcs = obj.dcs;
            model = obj.model;
            % Do Cross-Complementarity
            rbp = ~opts.right_boundary_point_explicit;

            % Impulse complementarities
            if opts.no_initial_impacts
                start_fe = 2;
            else
                start_fe = 1;
            end
            for ii=1:opts.N_stages
                for jj=start_fe:opts.N_finite_elements(ii)
                    Gij = {};
                    Hij = {};

                    Gij = [Gij,{obj.w.Lambda_normal(ii,jj)}];
                    Hij = [Hij,{obj.w.Y_gap(ii,jj) + obj.w.P_vn(ii,jj) + obj.w.N_vn(ii,jj)}];
                    Gij = [Gij,{obj.w.P_vn(ii,jj)}];
                    Hij = [Hij,{obj.w.N_vn(ii,jj)}];
                    if model.friction_exists
                        switch opts.friction_model
                          case 'Polyhedral'
                            Gij = [Gij,{obj.w.Delta_d(ii,jj)}];
                            Hij = [Hij,{obj.w.Lambda_tangent(ii,jj)}];
                            Gij = [Gij,{obj.w.Gamma_d(ii,jj)}];
                            Hij = [Hij,{obj.w.Beta_d(ii,jj)}];
                          case 'Conic'
                            Gij = [Gij,{obj.w.Gamma(ii,jj)}];
                            Hij = [Hij,{obj.w.Beta(ii,jj)}];
                            switch opts.conic_model_switch_handling
                              case 'Plain'
                              case 'Abs'
                                Gij = [Gij,{obj.w.P_vt(ii,jj)}];
                                Hij = [Hij,{obj.w.N_vt(ii,jj)}];
                              case 'Lp'
                                Gij = [Gij,{obj.w.Alpha_vt(ii,jj)}];
                                Hij = [Hij,{obj.w.P_vt(ii,jj)}];
                                Gij = [Gij,{1-obj.w.Alpha_vt(ii,jj)}];
                                Hij = [Hij,{obj.w.N_vt(ii,jj)}];
                            end
                        end
                    end
                    obj.G.impulse_comp(ii,jj) = {vertcat(Gij{:})};
                    obj.H.impulse_comp(ii,jj) = {vertcat(Hij{:})};
                end
            end
            % contact complementarities
            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};

                            % contact complementarities
                            if jj ~= 1 || ~opts.no_initial_impacts
                                for rr=1:opts.n_s
                                    lambda_normal_ijr = obj.w.lambda_normal(ii,jj,rr);

                                    Gij = vertcat(Gij, {lambda_normal_ijr});
                                    Hij = vertcat(Hij, {obj.w.Y_gap(ii,jj)});
                                end
                            end
                            for kk=1:(opts.n_s)
                                lambda_normal_ijk = obj.w.lambda_normal(ii,jj,kk);
                                for rr=1:(opts.n_s + rbp)
                                    y_gap_ijr = obj.w.y_gap(ii,jj,rr);
                                    
                                    Gij = vertcat(Gij, {lambda_normal_ijk});
                                    Hij = vertcat(Hij, {y_gap_ijr});
                                end
                            end
                            
                            if model.friction_exists
                                switch opts.friction_model
                                  case 'Polyhedral'
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        for rr=1:opts.n_s
                                            beta_d_ijr = obj.w.beta_d(ii,jj,rr);
                                            lambda_tangent_ijr = obj.w.lambda_tangent(ii,jj,rr);

                                            Gij = vertcat(Gij, {lambda_tangent_ijr});
                                            Hij = vertcat(Hij, {obj.w.Delta_d(ii,jj)});
                                            Gij = vertcat(Gij, {beta_d_ijr});
                                            Hij = vertcat(Hij, {obj.w.Gamma_d(ii,jj)});
                                        end
                                    end
                                    for kk=1:(opts.n_s)
                                        beta_d_ijk = obj.w.beta_d(ii,jj,kk);
                                        lambda_tangent_ijk = obj.w.lambda_tangent(ii,jj,kk);
                                        for rr=1:opts.n_s
                                            delta_d_ijr = obj.w.delta_d(ii,jj,rr);
                                            gamma_d_ijr = obj.w.gamma_d(ii,jj,rr);

                                            Gij = vertcat(Gij, {lambda_tangent_ijk});
                                            Hij = vertcat(Hij, {delta_d_ijr});
                                            Gij = vertcat(Gij, {beta_ijk});
                                            Hij = vertcat(Hij, {gamma_ijr});
                                        end
                                    end
                                  case 'Conic'
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        for rr=1:opts.n_s
                                            beta_ijr = obj.w.beta(ii,jj,rr);
                                            
                                            Gij = vertcat(Gij, {beta_ijr});
                                            Hij = vertcat(Hij, {obj.w.Gamma(ii,jj)});
                                        end
                                    end
                                    for kk=1:(opts.n_s)
                                        beta_ijk = obj.w.beta(ii,jj,kk);
                                        for rr=1:opts.n_s
                                            gamma_ijr = obj.w.gamma(ii,jj,rr);
                                            
                                            Gij = vertcat(Gij, {beta_ijk});
                                            Hij = vertcat(Hij, {gamma_ijr});
                                        end
                                    end
                                    switch opts.conic_model_switch_handling
                                      case 'Plain'
                                        % no extra expr
                                      case 'Abs'
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            for rr=1:opts.n_s
                                                p_vt_ijr = obj.w.p_vt(ii,jj,rr);
                                                
                                                Gij = vertcat(Gij, {p_vt_ijr});
                                                Hij = vertcat(Hij, {obj.w.N_vt(ii,jj)});
                                            end
                                            for rr=1:opts.n_s
                                                n_vt_ijr = obj.w.n_vt(ii,jj,rr);
                                                
                                                Gij = vertcat(Gij, {obj.w.P_vt(ii,jj)});
                                                Hij = vertcat(Hij, {n_vt_ijr});
                                            end
                                        end
                                        for kk=1:(opts.n_s)
                                            p_vt_ijk = obj.w.p_vt(ii,jj,kk);
                                            for rr=1:opts.n_s
                                                n_vt_ijr = obj.w.n_vt(ii,jj,rr);
                                                
                                                Gij = vertcat(Gij, {p_vt_ijk});
                                                Hij = vertcat(Hij, {n_vt_ijr});
                                            end
                                        end
                                      case 'Lp'
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            for rr=1:opts.n_s
                                                alpha_vt_ijr = obj.w.alpha_vt(ii,jj,rr);
                                                
                                                Gij = vertcat(Gij, {alpha_vt_ijr});
                                                Hij = vertcat(Hij, {obj.w.P_vt(ii,jj)});
                                                Gij = vertcat(Gij, {1-alpha_vt_ijr});
                                                Hij = vertcat(Hij, {obj.w.N_vt(ii,jj)});
                                            end
                                        end
                                        for kk=1:(opts.n_s)
                                            alpha_vt_ijk = obj.w.alpha_vt(ii,jj,kk);
                                            for rr=1:opts.n_s
                                                p_vt_ijr = obj.w.p_vt(ii,jj,rr);
                                                n_vt_ijr = obj.w.n_vt(ii,jj,rr);
                                                
                                                Gij = vertcat(Gij, {alpha_vt_ijk});
                                                Hij = vertcat(Hij, {p_vt_ijr});
                                                Gij = vertcat(Gij, {1-alpha_vt_ijk});
                                                Hij = vertcat(Hij, {n_vt_ijr});
                                            end
                                        end
                                    end
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                        end
                    end
                  case CrossCompMode.FE_STAGE
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            sum_lambda_normal = sum2(obj.w.lambda_normal(ii,jj,:));
                            if jj ~= 1 || ~opts.no_initial_impacts
                                Gij = vertcat(Gij, {sum_lambda_normal});
                                Hij = vertcat(Hij, {obj.w.Y_gap(ii,jj)});
                            end
                            for rr=1:(opts.n_s + rbp)
                                y_gap_ijr = obj.w.y_gap(ii,jj,rr);
                                Gij = vertcat(Gij, {sum_lambda_normal});
                                Hij = vertcat(Hij, {y_gap_ijr});
                            end
                            
                            if model.friction_exists
                                switch opts.friction_model
                                  case 'Polyhedral'
                                    sum_beta_d = sum2(obj.w.beta_d(ii,jj,:));
                                    sum_lambda_tangent = sum2(obj.w.lambda_tangent(ii,jj,:));
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        Gij = vertcat(Gij, {sum_lambda_tangent});
                                        Hij = vertcat(Hij, {obj.w.Delta_d(ii,jj)});
                                        Gij = vertcat(Gij, {sum_beta_d});
                                        Hij = vertcat(Hij, {obj.w.Gamma_d(ii,jj)}); 
                                    end
                                    for rr=1:opts.n_s
                                        delta_d_ijr = obj.w.delta_d(ii,jj,rr);
                                        gamma_d_ijr = obj.w.gamma_d(ii,jj,rr);

                                        Gij = vertcat(Gij, {sum_lambda_tangent});
                                        Hij = vertcat(Hij, {delta_d_ijr});
                                        Gij = vertcat(Gij, {sum_beta});
                                        Hij = vertcat(Hij, {gamma_ijr});
                                    end
                                  case 'Conic'
                                    sum_beta = sum2(obj.w.beta(ii,jj,:));
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        Gij = vertcat(Gij, {sum_beta});
                                        Hij = vertcat(Hij, {obj.w.Gamma(ii,jj)});
                                    end
                                    for rr=1:opts.n_s
                                        gamma_ijr = obj.w.gamma(ii,jj,rr);
                                        Gij = vertcat(Gij, {sum_beta});
                                        Hij = vertcat(Hij, {gamma_ijr});
                                    end
                                    switch opts.conic_model_switch_handling
                                      case 'Plain'
                                        % no extra expr
                                      case 'Abs'
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            p_vt_lb = obj.w.P_vt(ii,jj);
                                        else
                                            p_vt_lb = 0;
                                        end
                                        sum_p_vt = sum2(obj.w.p_vt(ii,jj,:)); 
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                                Gij = vertcat(Gij, {sum_p_vt});
                                                Hij = vertcat(Hij, {obj.w.N_vt(ii,jj)});
                                        end
                                        for rr=1:opts.n_s
                                            n_vt_ijr = obj.w.n_vt(ii,jj,rr);
                                            
                                            Gij = vertcat(Gij, {p_vt_lb + sum_p_vt});
                                            Hij = vertcat(Hij, {n_vt_ijr});
                                        end
                                      case 'Lp'
                                        sum_alpha_vt = sum2(obj.w.alpha_vt(ii,jj,:));
                                        sum_alpha_vt_minus = sum2(1-obj.w.alpha_vt(ii,jj,:));
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            Gij = vertcat(Gij, {sum_alpha_vt});
                                            Hij = vertcat(Hij, {obj.w.P_vt(ii,jj)});
                                            Gij = vertcat(Gij, {sum_alpha_vt_minus});
                                            Hij = vertcat(Hij, {obj.w.N_vt(ii,jj)});
                                        end
                                        for rr=1:opts.n_s
                                            p_vt_ijr = obj.w.p_vt(ii,jj,rr);
                                            n_vt_ijr = obj.w.n_vt(ii,jj,rr);
                                            
                                            Gij = vertcat(Gij, {sum_alpha_vt});
                                            Hij = vertcat(Hij, {p_vt_ijr});
                                            Gij = vertcat(Gij, {sum_alpha_vt_minus});
                                            Hij = vertcat(Hij, {n_vt_ijr});
                                        end
                                    end
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            sum_y_gap = obj.w.Y_gap(ii,jj) + sum2(obj.w.y_gap(ii,jj,:));
                            for kk=1:(opts.n_s + rbp)
                                lambda_normal_ijk = obj.w.lambda_normal(ii,jj,kk);
                                Gij = vertcat(Gij, {lambda_normal_ijk});
                                Hij = vertcat(Hij, {sum_y_gap});
                            end
                            if model.friction_exists
                                switch opts.friction_model
                                  case 'Polyhedral'
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        delta_d_lb = obj.w.Delta_d(ii,jj);
                                    else
                                        delta_d_lb = 0;
                                    end
                                    sum_delta_d = delta_d_lb + sum2(obj.w.delta_d(ii,jj,:));
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        gamma_d_lb = obj.w.Gamma_d(ii,jj);
                                    else
                                        gamma_d_lb = 0;
                                    end
                                    sum_gamma_d = gamma_d_lb + sum2(obj.w.gamma_d(ii,jj,:));
                                    
                                    for kk=1:opts.n_s
                                        lambda_tangent_ijk = obj.w.lambda_tangent(ii,jj,kk);
                                        beta_d_ijk = obj.w.beta_d(ii,jj,kk);

                                        Gij = vertcat(Gij, {lambda_tangent_ijk});
                                        Hij = vertcat(Hij, {sum_delta_d});
                                        Gij = vertcat(Gij, {beta_d_ijk});
                                        Hij = vertcat(Hij, {sum_gamma_d});
                                    end
                                  case 'Conic'
                                    sum_gamma = obj.w.Gamma(ii,jj) + sum2(obj.w.gamma(ii,jj,:));
                                    for kk=1:opts.n_s
                                        beta_ijk = obj.w.beta(ii,jj,kk);
                                        Gij = vertcat(Gij, {beta_ijk});
                                        Hij = vertcat(Hij, {sum_gamma});
                                    end
                                    switch opts.conic_model_switch_handling
                                      case 'Plain'
                                        % no extra expr
                                      case 'Abs'
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            n_vt_lb = obj.w.N_vt(ii,jj);
                                        else
                                            n_vt_lb = 0;
                                        end
                                        sum_n_vt = sum2(obj.w.n_vt(ii,jj,:));
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            Gij = vertcat(Gij, {obj.w.P_vt(ii,jj)});
                                            Hij = vertcat(Hij, {sum_n_vt});
                                        end
                                        for kk=1:opts.n_s
                                            p_vt_ijk = obj.w.n_vt(ii,jj,kk);
                                            
                                            Gij = vertcat(Gij, {p_vt_ijk});
                                            Hij = vertcat(Hij, {n_vt_lb + sum_n_vt});
                                        end
                                      case 'Lp'
                                        sum_p_vt = obj.w.P_vt(ii,jj) + sum2(obj.w.p_vt(ii,jj,:));
                                        sum_n_vt = obj.w.N_vt(ii,jj) + sum2(obj.w.n_vt(ii,jj,:));
                                        for kk=1:opts.n_s
                                            alpha_vt_ijk = obj.w.alpha_vt(ii,jj,kk);
                                            
                                            Gij = vertcat(Gij, {alpha_vt_ijk});
                                            Hij = vertcat(Hij, {sum_p_vt});
                                            Gij = vertcat(Gij, {1-alpha_vt_ijk});
                                            Hij = vertcat(Hij, {sum_n_vt});
                                        end
                                    end
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                        end
                    end
                  case CrossCompMode.FE_FE
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            if jj ~= 1 || ~opts.no_initial_impacts
                                y_gap_lb = obj.w.Y_gap(ii,jj);
                            else
                                y_gap_lb = 0;
                            end
                            sum_y_gap = y_gap_lb + sum2(obj.w.y_gap(ii,jj,:));
                            sum_lambda_normal = sum2(obj.w.lambda_normal(ii,jj,:));
                            Gij = vertcat(Gij, {sum_lambda_normal});
                            Hij = vertcat(Hij, {sum_y_gap});
                            if model.friction_exists
                                switch opts.friction_model
                                  case 'Polyhedral'
                                    sum_lambda_tangent = sum2(obj.w.lambda_tangent(ii,jj,:));
                                    sum_beta_d = sum2(obj.w.beta_d(ii,jj,:));
                                    
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        delta_d_lb = obj.w.Delta_d(ii,jj);
                                    else
                                        delta_d_lb = 0;
                                    end
                                    sum_delta_d = delta_d_lb + sum2(obj.w.delta_d(ii,jj,:));
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        gamma_d_lb = obj.w.Gamma_d(ii,jj);
                                    else
                                        gamma_d_lb = 0;
                                    end
                                    sum_gamma_d = gamma_d_lb + sum2(obj.w.gamma_d(ii,jj,:));

                                    Gij = vertcat(Gij, {sum_lambda_tangent});
                                    Hij = vertcat(Hij, {sum_delta_d});
                                    Gij = vertcat(Gij, {sum_beta_d});
                                    Hij = vertcat(Hij, {sum_gamma_d});
                                  case 'Conic'
                                    if jj ~= 1 || ~opts.no_initial_impacts
                                        gamma_lb = obj.w.Gamma(ii,jj);
                                    else
                                        gamma_lb = 0;
                                    end
                                    sum_gamma = gamma_lb + sum2(obj.w.gamma(ii,jj,:));
                                    sum_beta = sum2(obj.w.beta(ii,jj,:));
                                    Gij = vertcat(Gij, {sum_beta});
                                    Hij = vertcat(Hij, {sum_gamma});
                                    switch opts.conic_model_switch_handling
                                      case 'Plain'
                                        % no extra expr
                                      case 'Abs'
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            p_vt_lb = obj.w.P_vt(ii,jj);
                                        else
                                            p_vt_lb = 0;
                                        end
                                        sum_p_vt = sum2(obj.w.p_vt(ii,jj,:));
                                        if jj ~= 1 || ~opts.no_initial_impacts
                                            n_vt_lb = obj.w.N_vt(ii,jj);
                                        else
                                            n_vt_lb = 0;
                                        end
                                        sum_n_vt = sum2(obj.w.n_vt(ii,jj,:));
                                        Gij = vertcat(Gij, {p_vt_lb + sum_p_vt});
                                        Hij = vertcat(Hij, {n_vt_lb + sum_n_vt});
                                      case 'Lp'
                                        sum_p_vt = obj.w.P_vt(ii,jj) + sum2(obj.w.p_vt(ii,jj,:));
                                        sum_n_vt = obj.w.N_vt(ii,jj) + sum2(obj.w.n_vt(ii,jj,:));
                                        sum_alpha_vt = sum2(obj.w.alpha_vt(ii,jj,:));
                                        sum_alpha_vt_minus = sum2(1-obj.w.alpha_vt(ii,jj,:));
                                            
                                        Gij = vertcat(Gij, {sum_alpha_vt});
                                        Hij = vertcat(Hij, {sum_p_vt});
                                        Gij = vertcat(Gij, {sum_alpha_vt_minus});
                                        Hij = vertcat(Hij, {sum_n_vt});
                                    end
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                        end
                    end
                end
            else
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii);
                        for kk=1:opts.n_s
                            Gij = {};
                            Hij = {};

                            lambda_normal_ijk = obj.w.lambda_normal(ii,jj,kk);
                            y_gap_ijk = obj.w.y_gap(ii,jj,kk);
                            
                            Gij = vertcat(Gij, {lambda_normal_ijk});
                            Hij = vertcat(Hij, {y_gap_ijk});
                            
                            if model.friction_exists
                                switch opts.friction_model
                                  case 'Polyhedral'
                                    beta_d_ijk = obj.w.beta_d(ii,jj,kk);
                                    lambda_tangent_ijk = obj.w.lambda_tangent(ii,jj,kk);
                                    delta_d_ijk = obj.w.delta_d(ii,jj,kk);
                                    gamma_d_ijk = obj.w.gamma_d(ii,jj,kk);

                                    Gij = vertcat(Gij, {lambda_tangent_ijk});
                                    Hij = vertcat(Hij, {delta_d_ijk});
                                    Gij = vertcat(Gij, {beta_ijk});
                                    Hij = vertcat(Hij, {gamma_ijk});
                                    
                                  case 'Conic'
                                        beta_ijk = obj.w.beta(ii,jj,kk);
                                        gamma_ijr = obj.w.gamma(ii,jj,kk);
                                        
                                        Gij = vertcat(Gij, {beta_ijk});
                                        Hij = vertcat(Hij, {gamma_ijk});
                                    switch opts.conic_model_switch_handling
                                      case 'Plain'
                                        % no extra expr
                                      case 'Abs'
                                        p_vt_ijk = obj.w.p_vt(ii,jj,kk);
                                        n_vt_ijk = obj.w.n_vt(ii,jj,kk);
                                        
                                        Gij = vertcat(Gij, {p_vt_ijk});
                                        Hij = vertcat(Hij, {n_vt_ijk});
                                      case 'Lp'
                                        alpha_vt_ijk = obj.w.alpha_vt(ii,jj,kk);
                                        p_vt_ijk = obj.w.p_vt(ii,jj,kk);
                                        n_vt_ijk = obj.w.n_vt(ii,jj,kk);
                                        
                                        Gij = vertcat(Gij, {alpha_vt_ijk});
                                        Hij = vertcat(Hij, {p_vt_ijk});
                                        Gij = vertcat(Gij, {1-alpha_vt_ijk});
                                        Hij = vertcat(Hij, {n_vt_ijk});
                                    end
                                end
                            end
                            obj.G.standard_comp(ii,jj,kk) = {vertcat(Gij{:})};
                            obj.H.standard_comp(ii,jj,kk) = {vertcat(Hij{:})};
                        end
                    end
                end
            end
        end

        function generate_step_equilibration_constraints(obj)
            import casadi.*
            model = obj.model;
            opts = obj.opts;
            dims = obj.dcs.dims;
            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            if ~opts.use_fesd % do nothing
                return
            end
            rbp = ~opts.right_boundary_point_explicit;

            switch obj.opts.step_equilibration
              case StepEquilibrationMode.heuristic_mean
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii)
                        h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                        obj.f = obj.f + obj.p.rho_h_p()*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case StepEquilibrationMode.heuristic_diff
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements(ii)
                        h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                        obj.f = obj.f + obj.p.rho_h_p()*(obj.w.h(ii,jj)-obj.w.h(ii,jj-1))^2;
                    end
                end
              case StepEquilibrationMode.l2_relaxed_scaled
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            sigma_c_B = sigma_c_B + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);

                            sigma_c_F = sigma_c_F + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p() * tanh(eta/opts.step_equilibration_sigma) * delta_h.^2;
                    end
                                end
              case StepEquilibrationMode.l2_relaxed
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            sigma_c_B = sigma_c_B + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);

                            sigma_c_F = sigma_c_F + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p() * eta * delta_h.^2
                    end
                end
              case StepEquilibrationMode.direct
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            sigma_c_B = sigma_c_B + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);

                            sigma_c_F = sigma_c_F + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;

                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.g.step_equilibration(ii,jj) = {eta*delta_h, 0, 0};
                    end
                end
                %obj.eta_fun = Function('eta_fun', {obj.w.sym}, {eta_vec});
              case StepEquilibrationMode.direct_homotopy
                error("not currently implemented")
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            sigma_c_B = sigma_c_B + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);

                            sigma_c_F = sigma_c_F + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        homotopy_eq = [eta*delta_h - sigma;eta*delta_h + sigma];
                        obj.g.step_equilibration(ii,jj) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                end
                %obj.eta_fun = Function('eta_fun', {obj.w.sym}, {eta_vec});
              case StepEquilibrationMode.mlcp
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            sigma_c_B = sigma_c_B + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);

                            sigma_c_F = sigma_c_F + dcs.c_fun(x_ijk,z_ijk, v_global, p);
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        lambda_lambda = obj.w.lambda_lambda(ii,jj);
                        lambda_c = obj.w.lambda_c(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_F;pi_lambda-sigma_lambda_B;sigma_lambda_F+sigma_lambda_B-pi_lambda],0,inf};
                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_F;pi_c-sigma_c_B;sigma_c_F+sigma_c_B-pi_c],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_c-lambda_lambda;
                            B_max-pi_lambda;
                            B_max-pi_c];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1)],
                            [0*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.G.step_eq_kkt_max(ii,jj) = {[(B_max-pi_lambda);(B_max-pi_c)]};
                        obj.H.step_eq_kkt_max(ii,jj) = {[lambda_lambda;lambda_c]};
                        
                        % eta calculation
                        eta_const = [eta-pi_c;eta-pi_lambda;eta-pi_c-pi_lambda+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(dims.n_lambda,1);-inf*ones(dims.n_lambda,1);zeros(dims.n_lambda,1)],
                            [zeros(dims.n_lambda,1);zeros(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        M=obj.p.T()/opts.N_stages;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + (1/h0)*nu*M;
                            delta_h - (1/h0)*nu*M];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                    end
                end
            end
        end

        function populate_problem(obj)
            obj.create_variables();
            obj.generate_direct_transcription_constraints();
            obj.generate_complementarity_constraints();
            obj.generate_step_equilibration_constraints();

            obj.populated = true;
        end

        function create_solver(obj, solver_options, plugin)
            if ~obj.populated
                obj.populate_problem();
            end

            if ~exist('plugin')
                plugin = 'scholtes_ineq';
            end

            % Sort by indices to recover almost block-band structure.
            obj.w.sort_by_index();
            obj.g.sort_by_index();

            solver_options.assume_lower_bounds = true;

            obj.solver = nosnoc.solver.mpccsol('Mpcc solver', plugin, obj.to_casadi_struct(), solver_options);
        end

        function stats = solve(obj)
            opts = obj.opts;
            T_val = obj.p.T().val;

            for ii=1:opts.N_stages
                for jj=1:opts.N_finite_elements(ii)
                    % Recalculate ubh and lbh based on T_val
                    h0 = T_val/(opts.N_stages*opts.N_finite_elements(ii));
                    ubh = (1 + opts.gamma_h) * h0;
                    lbh = (1 - opts.gamma_h) * h0;
                    if opts.time_rescaling && ~opts.use_speed_of_time_variables
                        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                        ubh = (1+opts.gamma_h)*h0*opts.s_sot_max;
                        lbh = (1-opts.gamma_h)*h0/opts.s_sot_min;
                    end
                    obj.w.h(ii,jj).lb = lbh;
                    obj.w.h(ii,jj).ub = ubh;
                end
            end

            stats = solve@vdx.problems.Mpcc(obj);
        end

        function results = get_results_struct(obj)
            opts = obj.opts;
            model = obj.model;

            rbp = ~opts.right_boundary_point_explicit;
            
            if opts.right_boundary_point_explicit
                results.x = obj.discrete_time_problem.w.x(:,:,obj.opts.n_s).res;
                results.z = obj.discrete_time_problem.w.z(:,:,obj.opts.n_s).res;
                results.lambda = obj.discrete_time_problem.w.lambda(:,:,obj.opts.n_s).res;
                results.mu = obj.discrete_time_problem.w.mu(:,:,obj.opts.n_s).res;
                results.theta = obj.discrete_time_problem.w.theta(:,:,obj.opts.n_s).res;
            else
                results.x = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.z = [obj.discrete_time_problem.w.z(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.z(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.lambda = [obj.discrete_time_problem.w.lambda(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.lambda(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.mu = [obj.discrete_time_problem.w.mu(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.mu(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.theta = obj.discrete_time_problem.w.theta(:,:,obj.opts.n_s+1).res;
            end
            results.u = obj.w.u.res;

            % TODO also speed of time/etc here.
            if opts.use_fesd
                results.h = obj.w.h.res;
            else
                results.h = [];
                for ii=1:opts.N_stages
                    h = obj.p.T.val/(opts.N_stages*opts.N_finite_elements(ii));
                    results.h = [results.h,h*ones(1, opts.N_finite_elements(ii))];
                end
            end
            results.t_grid = cumsum([0, results.h]);
            if ~isempty(results.u)
                if opts.use_fesd
                    t_grid_u = [0];
                    for ii=1:opts.N_stages
                        h_sum = sum(obj.discrete_time_problem.w.h(ii,:).res);
                        t_grid_u = [t_grid_u, t_grid(end)+h_sum];
                    end
                    results.t_grid_u = t_grid_u;
                else
                    results.t_grid_u = linspace(0, obj.p.T.val, opts.N_stages+1);
                end
            end

            fields = fieldnames(S);
            empty_fields = cellfun(@(field) isempty(results.(field)), fields);
            results = rmfield(S, fields(empty_fields));
        end
    end
end
