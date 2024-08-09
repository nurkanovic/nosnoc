classdef Stewart < vdx.problems.Mpcc
    properties
        model % Model object
        dcs % Dcs object
        opts % Option object
        
        populated % Boolean variable indicating are the pupulating functions executed.
    end

    methods
        function obj = Stewart(dcs, opts)
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

            % OCP/Simulation problem parameters 
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
                    obj.p.rho_sot = {{'rho_sot',1}, opts.rho_sot};
                    obj.f = obj.f + obj.p.rho_sot()*sum((obj.w.sot(:)-1).^2);
                end
            end
            % Remark, VDX syntax for defining parameters
            % obj.p.parameter_name(indexing) = {{'parameter_name', length}, parameter_value};
            obj.p.p_time_var(1:opts.N_stages) = {{'p_time_var', dims.n_p_time_var}};
            for ii=1:opts.N_stages
                obj.p.p_time_var(ii).val = model.p_time_var_val(:,ii);
            end

            % TODO(@anton) This _severely_ hurts performance over the vectorized assignment by doing N_stages vertcats of
            %              casadi symbolics vs just a vectorized assignment which does one. As such there needs to be backend
            %              work done for vdx to cache vertcats of SX somehow. Current theory is one can simply keep a queue of
            %              symbolics to be added in a cell array until a read is done, at which point we call a single vertcat
            %              on the whole queue which is _significantly_ faster.
            % 2d vars: Variables that are defined for each finite element (and each control stage)
            if opts.use_fesd && opts.use_numerical_clock_state
                obj.w.numerical_time(0,0) = {{'t0', 1}};
                obj.w.numerical_time(0,0).lb = 0;
                obj.w.numerical_time(0,0).ub = 0;
            end
            for ii=1:opts.N_stages
                % other derived values
                h0 = opts.h_k(ii); % initial guess for length of current FE
                if obj.opts.use_fesd
                    ubh = (1 + opts.gamma_h) * h0; % upper bound for FE length
                    lbh = (1 - opts.gamma_h) * h0; % lower bound for FE length
                    if opts.time_rescaling && ~opts.use_speed_of_time_variables
                        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                        ubh = ubh*opts.s_sot_max;
                        lbh = lbh/opts.s_sot_min;
                    elseif opts.time_optimal_problem
                        ubh = ubh*(opts.T_final_max/opts.T);
                        lbh = lbh/((opts.T_final_min+eps)/opts.T);
                    end
                    % define finte elements lengths as variables
                    obj.w.h(ii,1:opts.N_finite_elements(ii)) = {{'h', 1}, lbh, ubh, h0};
                end
                if opts.use_fesd && opts.use_numerical_clock_state
                    obj.w.numerical_time(ii, 1:opts.N_finite_elements(ii)) = {{['t_' num2str(ii)], 1}};
                end
                
                if obj.opts.step_equilibration == StepEquilibrationMode.linear_complementarity
                    % define auxiliary variables needed to write step equlibration constraints
                    % Remark: only this step equilibration mode needs auxliary variables.
                    obj.w.B_max(ii,2:opts.N_finite_elements(ii)) = {{'B_max', dims.n_lambda},-inf,inf};
                    obj.w.pi_theta(ii,2:opts.N_finite_elements(ii)) = {{'pi_theta', dims.n_theta},-inf,inf};
                    obj.w.pi_lambda(ii,2:opts.N_finite_elements(ii)) = {{'pi_lambda', dims.n_lambda},-inf,inf};
                    obj.w.theta_mult(ii,2:opts.N_finite_elements(ii)) = {{'theta_mult', dims.n_theta},0,inf};
                    obj.w.lambda_mult(ii,2:opts.N_finite_elements(ii)) = {{'lambda_mult', dims.n_lambda},0,inf};
                    obj.w.eta(ii,2:opts.N_finite_elements(ii)) = {{'eta', dims.n_lambda},0,inf};
                    obj.w.nu(ii,2:opts.N_finite_elements(ii)) = {{'nu', 1},0,inf};
                end
            end

            % For c_n ~= 1 case
            rbp = ~opts.right_boundary_point_explicit;
            
            % 3d vars: Variables defined on each rk stage (and every finite element, and every control stage)
            
            % Remark on VDX syntax obj.w.x(index_of_control_stage,index_of_fe,index_of_rk_stage)
            % some are also defined at the initial point:
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}, model.x0, model.x0, model.x0};
            obj.w.z(0,0,opts.n_s) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
            obj.w.lambda(0,0,opts.n_s) = {{['lambda'], dims.n_lambda},0,inf, 1};
            obj.w.theta(0,0,opts.n_s) = {{['theta'], dims.n_theta},0, inf, 1/dims.n_theta};
            obj.w.mu(0,0,opts.n_s) = {{'mu', dims.n_mu},-inf,inf};
            % others are defined across all levels
            for ii=1:opts.N_stages
                if (opts.rk_representation == RKRepresentation.integral ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    % Remark on VDX syntax obj.w.x(ii,1:_NFE,1:n_s+rbp) - vectorized definition of variables
                    obj.w.x(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                else
                    obj.w.x(ii,1:opts.N_finite_elements(ii),opts.n_s+rbp) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                end
                if (opts.rk_representation == RKRepresentation.differential ||...
                    opts.rk_representation == RKRepresentation.differential_lift_x)
                    obj.w.v(ii,1:opts.N_finite_elements(ii),1:opts.n_s) = {{'v', dims.n_x}};
                end
                % TODO @Anton, at some point we might provide initial guesse for lambda,mu, theta
                obj.w.z(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
                obj.w.lambda(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'lambda', dims.n_lambda},0, inf, 1};
                obj.w.theta(ii,1:opts.N_finite_elements(ii),1:(opts.n_s)) = {{'theta', dims.n_theta},0, inf, 1/dims.n_theta};
                obj.w.mu(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'mu', dims.n_mu},-inf,inf};

                % Handle x_box settings (i.e., manage at which points the box constraints are enforced)
                if ~opts.x_box_at_stg && opts.rk_representation ~= RKRepresentation.differential
                    obj.w.x(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).ub = inf*ones(dims.n_x, 1);
                end

                if ~opts.x_box_at_fe
                    obj.w.x(ii,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(ii,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).ub = inf*ones(dims.n_x, 1);
                end
            end
        end
        % define core runge kutta equations
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
            lambda_0 = obj.w.lambda(0,0,opts.n_s);
            theta_0 = obj.w.theta(0,0,opts.n_s);
            mu_0 = obj.w.mu(0,0,opts.n_s);
            
            % Remark on VDX syntax, define symbolic expressions for
            % constraints g. the constraints are also grouped into, e.g. .z user algebraic functions), 
            % .algebraic (KKT conditions from from the Stewart DCS)
            obj.g.z(0,0,opts.n_s) = {dcs.g_z_fun(x_0, z_0, obj.w.u(1), v_global, [p_global;obj.p.p_time_var(1)])};
            obj.g.algebraic(0,0,opts.n_s) = {dcs.g_alg_fun(x_0, z_0, lambda_0, theta_0, mu_0, obj.w.u(1), v_global, [p_global;obj.p.p_time_var(1)])};
            
            x_prev = obj.w.x(0,0,opts.n_s); % last point of previous FE, needed for continuity conditions
            if opts.use_fesd && opts.use_numerical_clock_state
                t_prev = obj.w.numerical_time(0,0);
            end
            for ii=1:opts.N_stages
                h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(ii)); % TODO@Anton, why not h0 = opts.h_k(ii), easier to read.
                
                u_i = obj.w.u(ii); % read symbolc control variable of current control stage.
                p_stage = obj.p.p_time_var(ii); % read symbolc parameter of current control stage.
                p = [p_global;p_stage];
                if obj.opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                    s_sot = obj.w.sot(ii); % here, sot is a vector
                elseif obj.opts.use_speed_of_time_variables
                    s_sot = obj.w.sot(); % here, sot is a scalar
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
                    if opts.use_fesd && opts.use_numerical_clock_state
                        t_curr = obj.w.numerical_time(ii, jj);
                        if ~opts.time_freezing
                            obj.g.numerical_time_integrator(ii,jj) = {t_curr - (t_prev + s_sot*h)};
                        else
                            obj.g.numerical_time_integrator(ii,jj) = {t_curr - (t_prev + h)};
                        end
                        t_prev = t_curr;
                    end

                    switch opts.rk_representation
                      case RKRepresentation.integral
                        % In integral representation stage variables are states.
                        x_ij_end = opts.D_rk(1)*x_prev;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk); % k-th RK stage variable, in j-th finite element, in i-th control stage
                            z_ijk = obj.w.z(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            theta_ijk = obj.w.theta(ii,jj,kk);
                            mu_ijk = obj.w.mu(ii,jj,kk);

                            f_ijk = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);
                            q_ijk = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);
                            xk = opts.C_rk(1, kk+1) * x_prev;
                            for rr=1:opts.n_s
                                x_ijr = obj.w.x(ii,jj,rr);
                                xk = xk + opts.C_rk(rr+1, kk+1) * x_ijr;
                            end
                            % add stagewise constraint to vdx mpcc.
                            obj.g.dynamics(ii,jj,kk) = {h * f_ijk - xk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p)};

                            x_ij_end = x_ij_end + opts.D_rk(kk+1)*x_ijk;
                            
                            % add path constraints at every stage
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if ~opts.euler_cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.B_rk(kk+1)*h*q_ijk;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            lambda_ijk = obj.w.lambda(ii,jj,opts.n_s+1);
                            mu_ijk = obj.w.mu(ii,jj,opts.n_s+1);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_ijk, mu_ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
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
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            theta_ijk = obj.w.theta(ii,jj,kk);
                            mu_ijk = obj.w.mu(ii,jj,kk);

                            f_ijk = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);
                            q_ijk = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {f_ijk - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p)};
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if ~opts.euler_cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*q_ijk;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            lambda_ijk = obj.w.lambda(ii,jj,opts.n_s+1);
                            mu_ijk = obj.w.mu(ii,jj,opts.n_s+1);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_ijk, mu_ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                        else
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ij_end - obj.w.x(ii,jj,opts.n_s)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
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
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            theta_ijk = obj.w.theta(ii,jj,kk);
                            mu_ijk = obj.w.mu(ii,jj,kk);

                            f_ijk = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);
                            q_ijk = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {f_ijk - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                            obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p)};
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if ~opts.euler_cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*q_ijk;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            x_ijk = obj.w.x(ii,jj,opts.n_s+1);
                            z_ijk = obj.w.z(ii,jj,opts.n_s+1);
                            lambda_ijk = obj.w.lambda(ii,jj,opts.n_s+1);
                            mu_ijk = obj.w.mu(ii,jj,opts.n_s+1);

                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ijk - x_ij_end};
                            obj.g.lp_stationarity(ii,jj,opts.n_s+1) = {dcs.g_lp_stationarity_fun(x_ijk, z_ijk, lambda_ijk, mu_ijk, v_global, p)};
                            obj.g.z(ii,jj,opts.n_s+1) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                        else
                            obj.g.dynamics(ii,jj,opts.n_s+1) = {x_ij_end - obj.w.x(ii,jj,opts.n_s)};
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            % if path constraints are evaluated at control and FE grid points
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                        end
                    end
                    x_prev = obj.w.x(ii,jj,opts.n_s+rbp);
                end
                if ~opts.g_path_at_stg && ~opts.g_path_at_fe
                    % if path constraints are only evaluated at the control grid nodes
                    x_i = obj.w.x(ii, opts.N_finite_elements(ii), opts.n_s);
                    z_i = obj.w.z(ii, opts.N_finite_elements(ii), opts.n_s);
                    obj.g.path(ii) = {dcs.g_path_fun(x_i, z_i, u_i, v_global, p), model.lbg_path, model.ubg_path};
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
                if opts.euler_cost_integration
                    obj.f = obj.f + dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, u_i, v_global, p);
                end

                % Clock Constraints
                % handle numerical time
                if opts.use_fesd && opts.equidistant_control_grid
                    % Create a vdx relaxation struct which holds the information needed by vdx to automatically
                    % relax a constraint vector. It takes the type of relaxation, the name of the vdx.Variable
                    % to store the slacks, and the name of the vdx.Variable to store the relaxation parameter.
                    % This struct can be passed as the last element of the rhs of a vdx variable assignment.
                    % See vdx documentation on relaxation: TODO add when vdx docs updated.
                    relax_num_time_struct = vdx.RelaxationStruct(opts.relax_terminal_numerical_time.to_vdx, 's_numerical_time', 'rho_numerical_time');
                    if opts.time_optimal_problem && ~opts.time_freezing
                        % if time optimal and not time freezing we want T_final to be the rhs of the
                        % equidistant control grid constraints.
                        ecg_rhs = obj.w.T_final()/opts.N_stages;
                    else
                        % if time freezing or not time optimal we want T to be the rhs of the
                        % equidistant control grid constraints.
                        ecg_rhs = obj.p.T()/opts.N_stages;
                    end

                    if opts.use_numerical_clock_state
                        curr_t = obj.w.numerical_time(ii,opts.N_finite_elements(ii));
                        obj.g.equidistant_numerical_grid(ii) = {curr_t - ii*ecg_rhs, relax_num_time_struct};
                    else
                        if ~opts.time_freezing
                            obj.g.equidistant_numerical_grid(ii) = {s_sot*sum_h - ecg_rhs, relax_num_time_struct};
                        else
                            obj.g.equidistant_numerical_grid(ii) = {sum_h - ecg_rhs, relax_num_time_struct};
                        end
                    end
                    % vdx.RelaxationStruct.is_relaxed returns true if the relaxation is not `NONE`.
                    if relax_num_time_struct.is_relaxed
                        obj.p.rho_numerical_time().val = opts.rho_terminal_numerical_time;
                    end
                end

                % Handle possible physical time
                if opts.time_freezing && opts.stagewise_clock_constraint
                    x0 = obj.w.x(0,0,opts.n_s);
                    t0 = x0(end);
                    x_stage_end = obj.w.x(ii, opts.N_finite_elements(ii), opts.n_s+rbp);
                    t_stage_end = x_stage_end(end);

                    relax_phys_time_struct = vdx.RelaxationStruct(opts.relax_terminal_physical_time.to_vdx, 's_physical_time', 'rho_physical_time');
                    if opts.time_optimal_problem
                        obj.g.stagewise_clock_constraint(ii) = {t_stage_end - (ii*(obj.w.T_final()/opts.N_stages) + t0), relax_phys_time_struct};
                    else
                        obj.g.stagewise_clock_constraint(ii) = {t_stage_end - (ii*t_stage + t0), relax_phys_time_struct};
                    end
                    if relax_phys_time_struct.is_relaxed
                        obj.p.rho_physical_time().val = opts.rho_terminal_physical_time;
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
                x0 = obj.w.x(0,0,opts.n_s);
                t0 = x0(end);
                relax_phys_time_struct = vdx.RelaxationStruct(opts.relax_terminal_physical_time.to_vdx, 's_physical_time', 'rho_physical_time');
                % Terminal Phyisical Time (Possible terminal constraint on the clock state if time freezing is active).
                if opts.time_optimal_problem
                    % in the case of time optimal time freezing problem, physical clock state must reach final time.
                    obj.g.terminal_physical_time = {x_end(end)-(obj.w.T_final()+t0), relax_phys_time_struct};
                    if relax_phys_time_struct.is_relaxed
                        obj.p.rho_physical_time().val = opts.rho_terminal_physical_time;
                    end
                elseif opts.impose_terminal_phyisical_time && ~opts.stagewise_clock_constraint
                    % in the case of a generic time freezing problem, physical clock state must reach final time if
                    % we are imposing terminal physical time and we did not do so by apply stagewise constraints.
                    obj.g.terminal_physical_time = {x_end(end)-(obj.p.T()+t0), relax_phys_time_struct};
                    if relax_phys_time_struct.is_relaxed
                        obj.p.rho_physical_time().val = opts.rho_terminal_physical_time;
                    end
                end
                
            else
                if ~opts.use_fesd
                    if opts.time_optimal_problem
                        if opts.use_speed_of_time_variables
                            % without FESD we need speed of time variables and the sum of
                            % h*sot must be equal to the final time.
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
                            % not a time optimal problem so just sum of hs must be the terminal numerical time T.
                            relax_num_time_struct = vdx.RelaxationStruct(opts.relax_terminal_numerical_time.to_vdx, 's_numerical_time', 'rho_numerical_time');
                            obj.g.sum_h = {sum_h_all-obj.p.T(), relax_num_time_struct};
                            if relax_num_time_struct.is_relaxed
                                obj.p.rho_numerical_time().val = opts.rho_terminal_numerical_time;
                            end
                        else
                            if ~opts.use_speed_of_time_variables
                                % a time optimal problem without sot so sum of hs must be the optimal terminal numerical time T_final.
                                relax_num_time_struct = vdx.RelaxationStruct(opts.relax_terminal_numerical_time.to_vdx, 's_numerical_time', 'rho_numerical_time');
                                obj.g.sum_h = {sum_h_all-obj.w.T_final(), relax_num_time_struct};
                                if relax_num_time_struct.is_relaxed
                                    obj.p.rho_numerical_time().val = opts.rho_terminal_numerical_time;
                                end
                            else
                                % a time optimal problem with sot so sum of h*sot must be the optimal terminal "physical" time T_final.
                                % and the sum of h needs to be the terminal numerical time T
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

                                relax_num_time_struct = vdx.RelaxationStruct(opts.relax_terminal_numerical_time.to_vdx, 's_numerical_time', 'rho_numerical_time');
                                relax_phys_time_struct = vdx.RelaxationStruct(opts.relax_terminal_physical_time.to_vdx, 's_physical_time', 'rho_physical_time');
                                obj.g.sum_h = {sum_h_all-obj.p.T(), relax_num_time_struct};
                                obj.g.integral_clock = {sum_h_all-obj.w.T_final(), relax_phys_time_struct};
                                if relax_num_time_struct.is_relaxed
                                    obj.p.rho_numerical_time().val = opts.rho_terminal_numerical_time;
                                end
                                if relax_phys_time_struct.is_relaxed
                                    obj.p.rho_physical_time().val = opts.rho_terminal_physical_time;
                                end
                            end
                        end
                    end
                end
            end
            
            % Terminal constraint
            if opts.relax_terminal_constraint_homotopy
                error("nosnoc: Currently unsupported.")
            end
            g_terminal = dcs.g_terminal_fun(x_end, z_end, v_global, p_global);
            relax_terminal_struct = vdx.RelaxationStruct(opts.relax_terminal_constraint.to_vdx, 's_terminal', 'rho_terminal');
            obj.g.terminal = {g_terminal, model.lbg_terminal, model.ubg_terminal, relax_terminal_struct};
            if relax_terminal_struct.is_relaxed && all(size(g_terminal) > 0)
                obj.p.rho_terminal().val = opts.rho_terminal;
            end
        end

        function generate_complementarity_constraints(obj)
            import casadi.*
            opts = obj.opts;
            dcs = obj.dcs;
            model = obj.model;
            % Define standard complementarity and cross-complementarity constraints

            lambda_0 = obj.w.lambda(0,0,opts.n_s);
            theta_0 = obj.w.theta(0,0,opts.n_s);

            % inital_comp
            obj.G.initial_comp = {lambda_0};
            obj.H.initial_comp = {theta_0};
            
            rbp = ~opts.right_boundary_point_explicit;
            
            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            for rr=1:opts.n_s
                                theta_ijr = obj.w.theta(ii,jj,rr);

                                Gi = vertcat(Gij, {lambda_prev});
                                Hi = vertcat(Hij, {theta_ijr});
                            end    
                            for kk=1:(opts.n_s + rbp)
                                lambda_ijk = obj.w.lambda(ii,jj,kk);
                                for rr=1:opts.n_s
                                    theta_ijr = obj.w.theta(ii,jj,rr);

                                    Gij = vertcat(Gij, {lambda_ijk});
                                    Hij = vertcat(Hij, {theta_ijr});
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                        end
                    end
                  case CrossCompMode.FE_STAGE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            sum_theta = sum2(obj.w.theta(ii,jj,:));
                            Gij = {lambda_prev};
                            Hij = {sum_theta};
                            for kk=1:(opts.n_s + rbp)
                                lambda_ijk = obj.w.lambda(ii,jj,kk);

                                Gij = vertcat(Gij, {lambda_ijk});
                                Hij = vertcat(Hij, {sum_theta});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            sum_lambda = lambda_prev + sum2(obj.w.lambda(ii,jj,:));
                            Gij = {};
                            Hij = {};
                            for kk=1:opts.n_s
                                theta_ijk = obj.w.theta(ii,jj,kk);

                                Gij = vertcat(Gij, {sum_lambda});
                                Hij = vertcat(Hij, {theta_ijk});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                        end
                    end
                  case CrossCompMode.FE_FE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            sum_lambda = lambda_prev + sum2(obj.w.lambda(ii,jj,:));
                            sum_theta = sum2(obj.w.theta(ii,jj,:));
                            obj.G.cross_comp(ii,jj) = {sum_lambda};
                            obj.H.cross_comp(ii,jj) = {sum_theta};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                        end
                    end
                end
            else
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii);
                        Gij = {};
                        Hij = {};
                        for kk=1:opts.n_s
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            theta_ijk = obj.w.theta(ii,jj,kk);

                            obj.G.standard_comp(ii,jj, kk) = {lambda_ijk};
                            obj.H.standard_comp(ii,jj, kk) = {theta_ijk};
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
                        sigma_lambda_B = sum2(obj.w.lambda(ii,jj-1,:));
                        sigma_theta_B = sum2(obj.w.theta(ii,jj-1,:));
                        
                        sigma_lambda_F = obj.w.lambda(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda(ii,jj,:));
                        sigma_theta_F = sum2(obj.w.theta(ii,jj,:));
                        
                        pi_lambda = sigma_lambda_B .* sigma_lambda_F;
                        pi_theta = sigma_theta_B .* sigma_theta_F;
                        nu = pi_lambda + pi_theta;
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
                        sigma_lambda_B = sum2(obj.w.lambda(ii,jj-1,:));
                        sigma_theta_B = sum2(obj.w.theta(ii,jj-1,:));
                        
                        sigma_lambda_F = obj.w.lambda(ii,jj-1,opts.n_s + rbp) +sum2(obj.w.lambda(ii,jj,:));
                        sigma_theta_F = sum2(obj.w.theta(ii,jj,:));

                        pi_lambda = sigma_lambda_B .* sigma_lambda_F;
                        pi_theta = sigma_theta_B .* sigma_theta_F;
                        nu = pi_lambda + pi_theta;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p() * eta * delta_h.^2;
                    end
                end
              case StepEquilibrationMode.direct
                eta_vec = [];
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_B = sum2(obj.w.lambda(ii,jj-1,:));
                        sigma_theta_B = sum2(obj.w.theta(ii,jj-1,:));
                        
                        sigma_lambda_F = obj.w.lambda(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda(ii,jj,:));
                        sigma_theta_F = sum2(obj.w.theta(ii,jj,:));

                        pi_lambda = sigma_lambda_B .* sigma_lambda_F;
                        pi_theta = sigma_theta_B .* sigma_theta_F;
                        nu = pi_lambda + pi_theta;
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
                error("nosnoc: Currently not supported/implemented.")
                eta_vec = [];
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_lambda_B = sum2(obj.w.lambda(ii,jj-1,:));
                        sigma_theta_B = sum2(obj.w.theta(ii,jj-1,:));
                        
                        sigma_lambda_F = obj.w.lambda(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda(ii,jj,:));
                        sigma_theta_F = sum2(obj.w.theta(ii,jj,:));

                        pi_lambda = sigma_lambda_B .* sigma_lambda_F;
                        pi_theta = sigma_theta_B .* sigma_theta_F;
                        nu = pi_lambda + pi_theta;
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
                        sigma_lambda_B = sum2(obj.w.lambda(ii,jj-1,:));
                        sigma_theta_B = sum2(obj.w.theta(ii,jj-1,:));
                        
                        sigma_lambda_F = obj.w.lambda(ii,jj-1,opts.n_s + rbp) + sum2(obj.w.lambda(ii,jj,:));
                        sigma_theta_F = sum2(obj.w.theta(ii,jj,:));

                        lambda_mult = obj.w.lambda_mult(ii,jj);
                        theta_mult = obj.w.theta_mult(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        pi_theta = obj.w.pi_theta(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_F;pi_lambda-sigma_lambda_B;sigma_lambda_F+sigma_lambda_B-pi_lambda],0,inf};
                        obj.g.pi_theta_or(ii,jj) = {[pi_theta-sigma_theta_F;pi_theta-sigma_theta_B;sigma_theta_F+sigma_theta_B-pi_theta],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-theta_mult-lambda_mult;
                            B_max-pi_lambda;
                            B_max-pi_theta];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1)],
                            [0*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.G.step_eq_kkt_max(ii,jj) = {[(B_max-pi_lambda);(B_max-pi_theta)]};
                        obj.H.step_eq_kkt_max(ii,jj) = {[lambda_mult;theta_mult]};
                        
                        % eta calculation
                        eta_const = [eta-pi_theta;eta-pi_lambda;eta-pi_theta-pi_lambda+B_max];
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
                        elseif opts.time_optimal_problem
                            ubh = ubh*(opts.T_final_max/opts.T);
                            lbh = lbh/((opts.T_final_min+eps)/opts.T);
                        end
                        obj.w.h(ii,jj).lb = lbh;
                        obj.w.h(ii,jj).ub = ubh;
                    end
                end
            end

            stats = solve@vdx.problems.Mpcc(obj);
        end
    end
end
