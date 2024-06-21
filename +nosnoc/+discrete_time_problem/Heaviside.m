classdef Heaviside < vdx.problems.Mpcc
    properties
        model
        dcs
        opts

        populated
    end

    methods
        function obj = Heaviside(dcs, opts)
            obj = obj@vdx.problems.Mpcc();
            obj.model = dcs.model;
            obj.dcs = dcs;
            obj.opts = opts;
            obj.f = 0;
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
            % Remark: VDX syntax for defining variables: 
            % obj.w.variable = {{'variable_name', length}, lower_bound, upper_bound, initial_guess};
            obj.w.v_global = {{'v_global',dims.n_v_global}, model.lbv_global, model.ubv_global, model.v0_global};
            if opts.use_speed_of_time_variables && ~opts.local_speed_of_time_variable
                obj.w.sot = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
                if opts.time_freezing
                    obj.p.rho_sot = {{'rho_sot_p',1}, opts.rho_sot_p};
                    obj.f = obj.f + obj.p.rho_sot()*(obj.w.sot()-1)^2;
                end
            end
            if opts.time_optimal_problem
                obj.w.T_final = {{'T_final', 1}, opts.T_final_min, opts.T_final_max, opts.T};
                % obj.w.T_final(); gives back the symbolic variable of
                % T_final, the () indicates that it is a scalar CasADi symbolic variable
                obj.f = obj.f + obj.w.T_final();
            end

            % 1d vars: Variables that are defined per control stage
            obj.w.u(1:opts.N_stages) = {{'u', dims.n_u}, model.lbu, model.ubu, model.u0};
            if opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                obj.w.sot(1:opts.N_stages) = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
                if opts.time_freezing
                    obj.p.rho_sot = {{'rho_sot_p',1}, opts.rho_sot_p};
                    obj.f = obj.f + obj.p.rho_sot()*sum((obj.w.sot(:)-1).^2);
                end
            end
            % Remark, VDX syntax for defining parameters
            % obj.p.parameter_name(indexing) = {{'parameter_name', length}, parameter_value};
            obj.p.p_time_var(1:opts.N_stages) = {{'p_time_var', dims.n_p_time_var}, model.p_time_var_val};

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
                if obj.opts.step_equilibration == StepEquilibrationMode.linear_complementarity
                    obj.w.B_max(ii,2:opts.N_finite_elements(ii)) = {{'B_max', dims.n_lambda},-inf,inf};
                    obj.w.pi_lambda_n(ii,2:opts.N_finite_elements(ii)) = {{'pi_lambda_n', dims.n_lambda},-inf,inf};
                    obj.w.pi_lambda_p(ii,2:opts.N_finite_elements(ii)) = {{'pi_lambda_p', dims.n_lambda},-inf,inf};
                    obj.w.lambda_lambda_n(ii,2:opts.N_finite_elements(ii)) = {{'lambda_lambda_n', dims.n_lambda},0,inf};
                    obj.w.lambda_lambda_p(ii,2:opts.N_finite_elements(ii)) = {{'lambda_lambda_p', dims.n_lambda},0,inf};
                    obj.w.eta(ii,2:opts.N_finite_elements(ii)) = {{'eta', dims.n_lambda},0,inf};
                    obj.w.nu(ii,2:opts.N_finite_elements(ii)) = {{'nu', 1},0,inf};
                end
            end

            % For c_n ~= 1 case
            rbp = ~opts.right_boundary_point_explicit;
            
            % 3d vars: Variables defined on each rk stage
            %          some of which are also defined at the initial point:
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}, model.x0, model.x0, model.x0};
            obj.w.z(0,0,opts.n_s) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
            obj.w.lambda_n(0,0,opts.n_s) = {{['lambda_n'], dims.n_lambda},0,inf};
            obj.w.lambda_p(0,0,opts.n_s) = {{['lambda_p'], dims.n_lambda},0,inf};
            for ii=1:opts.N_stages
                if (opts.rk_representation == RKRepresentation.integral ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    obj.w.x(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                else
                    obj.w.x(ii,1:opts.N_finite_elements(ii),opts.n_s) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                end
                if (opts.rk_representation == RKRepresentation.differential ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    obj.w.v(ii,1:opts.N_finite_elements(ii),1:opts.n_s) = {{'v', dims.n_x}};
                end
                
                obj.w.z(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
                obj.w.alpha(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'alpha', dims.n_alpha},0, 1};
                obj.w.lambda_n(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'lambda_n', dims.n_lambda},0, inf};
                obj.w.lambda_p(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'lambda_p', dims.n_lambda},0, inf};
                
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
            lambda_n_0 = obj.w.lambda_n(0,0,opts.n_s);
            lambda_p_0 = obj.w.lambda_p(0,0,opts.n_s);

            obj.g.z(0,0,opts.n_s) = {dcs.g_z_fun(x_0, z_0, obj.w.u(1), v_global, [p_global;obj.p.p_time_var(1)])};
            obj.g.lp_stationarity(0,0,opts.n_s) = {dcs.g_lp_stationarity_fun(x_0, z_0, lambda_n_0, lambda_p_0, v_global, [p_global;obj.p.p_time_var(1)])};
            
            x_prev = obj.w.x(0,0,opts.n_s);
            for ii=1:opts.N_stages
                h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(ii));
                
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
                if opts.time_optimal_problem && ~opts.use_speed_of_time_variables
                    t_stage = obj.w.T_final()/(opts.N_stages*opts.N_finite_elements(ii));
                elseif opts.time_optimal_problem
                    t_stage = s_sot*obj.p.T()/opts.N_stages;
                else
                    t_stage = obj.p.T()/opts.N_stages;
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
                    
                    switch opts.rk_representation
                      case RKRepresentation.integral
                        % In integral representation stage variables are states.
                        x_ij_end = x_prev;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            alpha_ijk = obj.w.alpha(ii,jj,kk);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);
                            
                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);
                            xk = opts.C_rk(1, kk+1) * x_prev;
                            for rr=1:opts.n_s
                                x_ijr = obj.w.x(ii,jj,rr);
                                xk = xk + opts.C_rk(rr+1, kk+1) * x_ijr;
                            end
                            obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p)};

                            x_ij_end = x_ij_end + opts.D_rk(kk+1)*x_ijk;
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.B_rk(kk+1)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,opts.n_s+1);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,opts.n_s+1);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_p_ijk, lambda_p_ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                        end
                      case RKRepresentation.differential
                        % In differential representation stage variables are the state derivatives.
                        X_ijk = {};
                        for kk = 1:opts.n_s
                            x_temp = x_prev;
                            for rr = 1:opts.n_s
                                x_temp = x_temp + h*opts.A_rk(kk,rr)*obj.w.v(ii,jj,rr);
                            end
                            X_ijk = [X_ijk {x_temp}];
                        end
                        X_ijk = [X_ijk, {obj.w.x(ii,jj,opts.n_s)}];
                        x_ij_end = x_prev;
                        for kk=1:opts.n_s
                            x_ijk = X_ijk{kk};
                            v_ijk = obj.w.v(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            alpha_ijk = obj.w.alpha(ii,jj,kk);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);

                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {fj - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p)};
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            alpha_ijk = obj.w.alpha(ii,jj,kk);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_n_ijk, lambda_p_ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
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
                            x_temp = x_prev;
                            for rr = 1:opts.n_s
                                x_temp = x_temp + h*opts.A_rk(kk,rr)*obj.w.v(ii,jj,rr);
                            end
                            obj.g.lift_x(ii,jj,kk) = {x_ijk - x_temp};
                        end
                        x_ij_end = x_prev;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            v_ijk = obj.w.v(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            alpha_ijk = obj.w.alpha(ii,jj,kk);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);

                            fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);
                            qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {fj - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_p_ijk, ui, v_global, p)};
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, ui, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if opts.cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*qj;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_n_ijk, lambda_p.ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
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
                if ~isempty(model.x_ref_val)
                    obj.f = obj.f + t_stage*dcs.f_lsq_x_fun(obj.w.x(ii,opts.N_finite_elements(ii),opts.n_s),...
                        model.x_ref_val(:,ii),...
                        p);
                end
                if ~isempty(model.u_ref_val)
                    obj.f = obj.f + t_stage*dcs.f_lsq_u_fun(obj.w.u(ii),...
                        model.u_ref_val(:,ii),...
                        p);
                end
                if ~opts.cost_integration
                    obj.f = obj.f + dcs.f_q_fun(x_ijk, z_ijk, alpha_ijk, lambda_n_ijk, lambda_n_ijk, ui, v_global, p);
                end

                % Clock <Constraints
                % TODO(@anton) HERE BE DRAGONS. This is by far the worst part of current nosnoc as it requires the discrete problem
                %              to understand something about the time-freezing reformulation which is ugly.
                if opts.use_fesd && opts.equidistant_control_grid
                    if opts.time_optimal_problem
                        if opts.use_speed_of_time_variables
                            obj.g.equidistant_control_grid(ii) = {[sum_h - opts.h;s_sot*sum_h - obj.w.T_final()/opts.N_stages]};
                        else
                            obj.g.equidistant_control_grid(ii) = {sum_h - obj.w.T_final()/opts.N_stages};
                        end
                    elseif ~opts.time_freezing
                        if opts.relax_terminal_numerical_time
                            % Negative values sometimes help convergence.
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

                if opts.time_freezing && opts.stagewise_clock_constraint
                    x0 = obj.w.x(0,0,opts.n_s);
                    t0 = x0(end);
                    x_stage_end = obj.w.x(ii, opts.N_finite_elements(ii), opts.n_s+rbp);
                    t_stage_end = x_stage_end(end);
                    if problem_options.time_optimal_problem
                        obj.g.stagewise_clock_constraint(ii) = {t_stage_end - ii*(obj.w.T_final()/problem_options.N_stages) + t0};
                    else
                        obj.g.stagewise_clock_constraint(ii) = {t_stage_end - ii*t_stage + t0};
                    end
                end
            end

            x_end = obj.w.x(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            z_end = obj.w.z(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            % Terminal cost
            obj.f = obj.f + dcs.f_q_T_fun(x_end, z_end, v_global, p_global);
            
            % Terminal_lsq_cost
            if ~isempty(model.x_ref_end_val)
                obj.f = obj.f + h0*opts.N_finite_elements(ii)*dcs.f_lsq_T_fun(x_end,...
                    model.x_ref_end_val,...
                    p_global);
            end

            %  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
            % If the control grid is not equidistant, the constraint on sum of h happen only at the end.
            % The constraints are split to those which are related to numerical and physical time, to make it easier to read.
            if opts.time_freezing
                % Terminal Phyisical Time (Possible terminal constraint on the clock state if time freezing is active).
                if opts.time_optimal_problem
                    obj.g.terminal_physical_time = {x_end(end)-obj.w.T_final()};
                else
                    if opts.impose_terminal_phyisical_time && ~opts.stagewise_clock_constraint
                        obj.g.terminal_physical_time = {x_end(end)-obj.p.T()};
                    else
                        % no terminal constraint on the numerical time
                    end
                end
            else
                if ~opts.use_fesd
                    if opts.time_optimal_problem
                        % if time_freezing is on, everything is done via the clock state.
                        if opts.use_speed_of_time_variables
                            integral_clock_state = 0;
                            for ii=1:opts.N_stages
                                if opts.local_speed_of_time_variable
                                    s_sot = obj.d.sot(ii);
                                else
                                    s_sot = obj.w.sot();
                                end
                                for jj=1:opts.N_finite_elements(ii)
                                    h = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                                    integral_clock_state = integral_clock_state + h*s_sot;
                                end
                            end
                            obj.g.integral_clock_state = {integral_clock_state-obj.w.T_final()};
                        else
                            % otherwise treated via variable h_ki, i.e.,  h_ki =  T_final/(N_stages*N_FE)
                        end
                    end
                else
                    % if equidistant_control_grid = true all time constraint are added in
                    % the main control loop for every control stage k and the code
                    % below is skipped
                    if  ~opts.equidistant_control_grid
                        % T_num = T_phy = T_final =  T.
                        % all step sizes add up to prescribed time T.
                        % if use_speed_of_time_variables = true, numerical time is decupled from the sot scaling (no mather if local or not):
                        sum_h_all = sum2(obj.w.h(:,:));
                        
                        if ~opts.time_optimal_problem
                            if opts.relax_terminal_numerical_time
                                obj.w.s_numerical_time = {{'s_numerical', 1}, -2*opts.T, 2*opts.T, opts.T/2};
                                s_numerical = obj.w.s_numerical_time();
                                g_sum_h = [sum_h_all-obj.p.T()-s_numerical;-(sum_h_all-obj.p.T())-s_numerical]
                                obj.g.sum_h = {g_sum_h,-inf,0};
                                if opts.relax_terminal_numerical_time_homotopy
                                    %rho_terminal_numerical_time = 1/sigma_p;
                                    %obj.augmented_objective = obj.augmented_objective + rho_terminal_numerical_time*s_numerical;
                                    error('homotopy parameter not acessible here.')
                                else
                                    % TODO opts.rho_terminal_numerical_time as param
                                    obj.f = obj.f + opts.rho_terminal_numerical_time*s_numerical;
                                end
                            else
                                obj.g.sum_h = {sum_h_all-obj.p.T()};
                            end
                        else
                            if ~opts.use_speed_of_time_variables
                                if opts.relax_terminal_numerical_time
                                    obj.w.s_numerical_time = {{'s_numerical', 1}, -2*opts.T, 2*opts.T, opts.T/2};
                                    s_numerical = obj.w.s_numerical_time();
                                    g_sum_h = [sum_h_all-obj.w.T_final()-s_numerical;-(sum_h_all-obj.w.T_final())-s_numerical]
                                    obj.g.sum_h = {g_sum_h,-inf,0};
                                    if opts.relax_terminal_numerical_time_homotopy
                                        %rho_terminal_numerical_time = 1/sigma_p;
                                        %obj.augmented_objective = obj.augmented_objective + rho_terminal_numerical_time*s_numerical;
                                        error('homotopy parameter not acessible here.')
                                    else
                                        % TODO opts.rho_terminal_numerical_time as param
                                        obj.f = obj.f + opts.rho_terminal_numerical_time*s_numerical;
                                    end
                                else
                                    obj.g.sum_h = {sum_h_all-obj.w.T_final()};
                                end
                            else
                                integral_clock_state = 0;
                                for ii=1:opts.N_stages
                                    if opts.local_speed_of_time_variable
                                        s_sot = obj.d.sot(ii);
                                    else
                                        s_sot = obj.w.sot();
                                    end
                                    for jj=1:opts.N_finite_elements(ii)
                                        h = obj.w.h(ii,jj);
                                        integral_clock_state = integral_clock_state + h*s_sot;
                                    end
                                end
                                % T_num = T_phy = T_final \neq T.
                                if opts.relax_terminal_numerical_time
                                    obj.w.s_numerical_time = {{'s_numerical', 1}, -2*opts.T, 2*opts.T, opts.T/2};
                                    s_numerical = obj.w.s_numerical_time();
                                    g_integral_clock = [integral_clock_state-obj.w.T_final()-s_numerical;-(integral_clock_state-obj.w.T_final())-s_numerical]
                                    obj.g.integral_clock = {g_sum_h,-inf,0};
                                    if opts.relax_terminal_numerical_time_homotopy
                                        %rho_terminal_numerical_time = 1/sigma_p;
                                        %obj.augmented_objective = obj.augmented_objective + rho_terminal_numerical_time*s_numerical;
                                        error('homotopy parameter not acessible here.')
                                    else
                                        % TODO opts.rho_terminal_numerical_time as param
                                        obj.f = obj.f + opts.rho_terminal_numerical_time*s_numerical;
                                    end
                                else
                                    obj.g.sum_h = {sum_h_all-obj.p.T()};
                                    obj.g.integral_clock = {sum_h_all-obj.w.T_final()};
                                end
                            end
                        end
                    end
                end
            end
            
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

        function generate_complementarity_constraints(obj)
            import casadi.*
            opts = obj.opts;
            dcs = obj.dcs;
            model = obj.model;
            % Do Cross-Complementarity

            rbp = ~opts.right_boundary_point_explicit;
            
            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    lambda_n_prev = obj.w.lambda_n(0,0,opts.n_s);
                    lambda_p_prev = obj.w.lambda_p(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            for rr=1:opts.n_s
                                alpha_ijr = obj.w.alpha(ii,jj,rr);

                                Gij = vertcat(Gij, {lambda_n_prev}, {lambda_p_prev});
                                Hij = vertcat(Hij, {alpha_ijr}, {1-alpha_ijr});
                            end
                            for kk=1:(opts.n_s + rbp)
                                lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                                lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);
                                for rr=1:opts.n_s
                                    alpha_ijr = obj.w.alpha(ii,jj,rr);

                                    Gij = vertcat(Gij, {lambda_n_ijk}, {lambda_p_ijk});
                                    Hij = vertcat(Hij, {alpha_ijr}, {1-alpha_ijr});
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_n_prev = obj.w.lambda_n(ii,jj,opts.n_s+rbp);
                            lambda_p_prev = obj.w.lambda_p(ii,jj,opts.n_s+rbp);
                        end
                    end
                  case CrossCompMode.FE_STAGE                    
                    lambda_n_prev = obj.w.lambda_n(0,0,opts.n_s);
                    lambda_p_prev = obj.w.lambda_p(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii)
                            sum_alpha = sum2(obj.w.alpha(ii,jj,:));
                            sum_alpha_n = sum2(1-obj.w.alpha(ii,jj,:));
                            Gij = {lambda_n_prev; lambda_p_prev};
                            Hij = {sum_alpha; sum_alpha_n};
                            for kk=1:(opts.n_s + rbp)
                                lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                                lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);

                                Gij = vertcat(Gij, {lambda_n_ijk}, {lambda_p_ijk});
                                Hij = vertcat(Hij, {sum_alpha}, {sum_alpha_n});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_n_prev = obj.w.lambda_n(ii,jj,opts.n_s+rbp);
                            lambda_p_prev = obj.w.lambda_p(ii,jj,opts.n_s+rbp);
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    lambda_n_prev = obj.w.lambda_n(0,0,opts.n_s);
                    lambda_p_prev = obj.w.lambda_p(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii)
                            sum_lambda_n = lambda_n_prev + sum2(obj.w.lambda_n(ii,jj,:));
                            sum_lambda_p = lambda_p_prev + sum2(obj.w.lambda_p(ii,jj,:));
                            Gij = {};
                            Hij = {};
                            for kk=1:opts.n_s
                                alpha_ijk = obj.w.alpha(ii,jj,kk);

                                Gij = vertcat(Gij, {sum_lambda_n}, {sum_lambda_p});
                                Hij = vertcat(Hij, {alpha_ijk}, {1-alpha_ijk});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_n_prev = obj.w.lambda_n(ii,jj,opts.n_s+rbp);
                            lambda_p_prev = obj.w.lambda_p(ii,jj,opts.n_s+rbp);
                        end
                    end
                  case CrossCompMode.FE_FE
                    lambda_n_prev = obj.w.lambda_n(0,0,opts.n_s);
                    lambda_p_prev = obj.w.lambda_p(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii)
                            sum_lambda_n = lambda_n_prev + sum2(obj.w.lambda_n(ii,jj,:));
                            sum_lambda_p = lambda_p_prev + sum2(obj.w.lambda_p(ii,jj,:));
                            sum_alpha = sum2(obj.w.alpha(ii,jj,:));
                            sum_alpha_n = sum2(1-obj.w.alpha(ii,jj,:));
                            obj.G.cross_comp(ii,jj) = {[sum_lambda_n;sum_lambda_p]};
                            obj.H.cross_comp(ii,jj) = {[sum_alpha;sum_alpha_n]};
                            lambda_n_prev = obj.w.lambda_n(ii,jj,opts.n_s+rbp);
                            lambda_p_prev = obj.w.lambda_p(ii,jj,opts.n_s+rbp);
                        end
                    end
                end
            else
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii)
                        Gij = {};
                        Hij = {};
                        for kk=1:opts.n_s
                            lambda_n_ijk = obj.w.lambda_n(ii,jj,kk);
                            lambda_p_ijk = obj.w.lambda_p(ii,jj,kk);
                            alpha_ijk = obj.w.alpha(ii,jj,kk);

                            obj.G.standard_comp(ii,jj, kk) = {[lambda_n_ijk;lambda_p_ijk]};
                            obj.H.standard_comp(ii,jj, kk) = {[alpha_ijk;1-alpha_ijk]};
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
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_n_B = sum2(obj.w.lambda_n(ii,jj-1,:));
                        sigma_lambda_p_B = sum2(obj.w.lambda_p(ii,jj-1,:));
                        
                        sigma_lambda_n_F = obj.w.lambda_n(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_n(ii,jj,:));
                        sigma_lambda_p_F = obj.w.lambda_p(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_p(ii,jj,:));
                        
                        pi_lambda_n = sigma_lambda_n_B .* sigma_lambda_n_F;
                        pi_lambda_p = sigma_lambda_p_B .* sigma_lambda_p_F;
                        nu = pi_lambda_n + pi_lambda_p;
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
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_n_B = sum2(obj.w.lambda_n(ii,jj-1,:));
                        sigma_lambda_p_B = sum2(obj.w.lambda_p(ii,jj-1,:));
                        
                        sigma_lambda_n_F = obj.w.lambda_n(ii,jj-1,opts.n_s + rbp) +sum2(obj.w.lambda_n(ii,jj,:));
                        sigma_lambda_p_F = obj.w.lambda_p(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_p(ii,jj,:));

                        pi_lambda_n = sigma_lambda_n_B .* sigma_lambda_n_F;
                        pi_lambda_p = sigma_lambda_p_B .* sigma_lambda_p_F;
                        nu = pi_lambda_n + pi_lambda_p;
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
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_n_B = sum2(obj.w.lambda_n(ii,jj-1,:));
                        sigma_lambda_p_B = sum2(obj.w.lambda_p(ii,jj-1,:));
                        
                        sigma_lambda_n_F = obj.w.lambda_n(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_n(ii,jj,:));
                        sigma_lambda_p_F = obj.w.lambda_p(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_p(ii,jj,:));

                        pi_lambda_n = sigma_lambda_n_B .* sigma_lambda_n_F;
                        pi_lambda_p = sigma_lambda_p_B .* sigma_lambda_p_F;
                        nu = pi_lambda_n + pi_lambda_p;
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
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_n_B = sum2(obj.w.lambda_n(ii,jj-1,:));
                        sigma_lambda_p_B = sum2(obj.w.lambda_p(ii,jj-1,:));
                        
                        sigma_lambda_n_F = obj.w.lambda_n(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_n(ii,jj,:));
                        sigma_lambda_p_F = obj.w.lambda_p(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_p(ii,jj,:));

                        pi_lambda_n = sigma_lambda_n_B .* sigma_lambda_n_F;
                        pi_lambda_p = sigma_lambda_p_B .* sigma_lambda_p_F;
                        nu = pi_lambda_n + pi_lambda_p;
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
              case StepEquilibrationMode.linear_complementarity
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements(ii)
                        h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                        sigma_lambda_n_B = sum2(obj.w.lambda_n(ii,jj-1,:));
                        sigma_lambda_p_B = sum2(obj.w.lambda_p(ii,jj-1,:));
                        
                        sigma_lambda_n_F = obj.w.lambda_n(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_n(ii,jj,:));
                        sigma_lambda_p_F = obj.w.lambda_p(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda_p(ii,jj,:));

                        lambda_lambda_n = obj.w.lambda_lambda_n(ii,jj);
                        lambda_lambda_p = obj.w.lambda_lambda_p(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_lambda_n = obj.w.pi_lambda_n(ii,jj);
                        pi_lambda_p = obj.w.pi_lambda_p(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_lambda_n_or(ii,jj) = {[pi_lambda_n-sigma_lambda_n_F;pi_lambda_n-sigma_lambda_n_B;sigma_lambda_n_F+sigma_lambda_n_B-pi_lambda_n],0,inf};
                        obj.g.pi_lambda_p_or(ii,jj) = {[pi_lambda_p-sigma_lambda_p_F;pi_lambda_p-sigma_lambda_p_B;sigma_lambda_p_F+sigma_lambda_p_B-pi_lambda_p],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_lambda_n-lambda_lambda_p;
                            B_max-pi_lambda_n;
                            B_max-pi_lambda_p];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1)],
                            [0*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.G.step_eq_kkt_max(ii,jj) = {[(B_max-pi_lambda_n);(B_max-pi_lambda_p)]};
                        obj.H.step_eq_kkt_max(ii,jj) = {[lambda_lambda_n;lambda_lambda_p]};
                        
                        % eta calculation
                        eta_const = [eta-pi_lambda_p;eta-pi_lambda_n;eta-pi_lambda_p-pi_lambda_n+B_max];
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

            if opts.use_fesd
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
                results.lambda_n = obj.discrete_time_problem.w.lambda_n(:,:,obj.opts.n_s).res;
                results.lambda_p = obj.discrete_time_problem.w.lambda_p(:,:,obj.opts.n_s).res;
                results.alpha = obj.discrete_time_problem.w.alpha(:,:,obj.opts.n_s).res;
            else
                results.x = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.z = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.lambda_n = [obj.discrete_time_problem.w.lambda_n(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.lambda_n(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.lambda_p = [obj.discrete_time_problem.w.lambda_p(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.lambda_p(1:opts.N_stages,:,obj.opts.n_s+1).res];
                results.alpha = obj.discrete_time_problem.w.alpha(:,:,obj.opts.n_s+1).res;
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
