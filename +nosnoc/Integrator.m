classdef Integrator < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_time_problem
        stats

        solver_exists = false;
    end

    properties(Access=private)
        w_all
    end

    methods
        function obj = Integrator(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            % Always process model and options
            % for integrator also take the extra step of re-calculating N_stages/N_finite_elements
            opts.preprocess();
            if opts.N_stages > 1
                warning("Integrator created with more than 1 control stage. Converting this to finite elements.")
                N_fe = sum(opts.N_finite_elements);
                opts.N_finite_elements = N_fe;
                opts.N_stages = 1;
            end
            model.verify_and_backfill(opts);

            if class(model) == "nosnoc.model.Cls" && opts.time_freezing
                model = nosnoc.time_freezing.reformulation(model, opts);
                model.verify_and_backfill(opts);
                obj.model = model;
            end

            % Run pipeline
            switch class(model)
              case "nosnoc.model.Pss"
                if opts.dcs_mode == DcsMode.Stewart
                    obj.dcs = nosnoc.dcs.Stewart(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_time_problem.Stewart(obj.dcs, opts);
                    obj.discrete_time_problem.populate_problem();
                elseif opts.dcs_mode == DcsMode.Heaviside
                    obj.dcs = nosnoc.dcs.Heaviside(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_time_problem.Heaviside(obj.dcs, opts);
                    obj.discrete_time_problem.populate_problem();
                else
                    error("PSS models can only be reformulated using the Stewart or Heaviside Step reformulations.")
                end
              case "nosnoc.model.Heaviside"
                obj.dcs = nosnoc.dcs.Heaviside(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Heaviside(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Cls"
                obj.dcs = nosnoc.dcs.Cls(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Cls(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Pds"
                if ~opts.right_boundary_point_explicit
                    error("You are using an rk scheme with its right boundary point (c_n) not equal to one. Please choose another scheme e.g. RADAU_IIA")
                end
                if opts.rk_representation == RKRepresentation.differential
                    error("Differential representation without lifting is unsupported for gradient complementarity systems. Use integral or lifted differential representation")
                end
                obj.dcs = nosnoc.dcs.Gcs(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Gcs(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
            end
        end

        function stats = solve(obj, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            if ~obj.solver_exists
                obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);
                if isfield(obj.solver_opts.opts_casadi_nlp,'iteration_callback')
                    obj.solver_opts.opts_casadi_nlp.iteration_callback.mpcc = obj.discrete_time_problem;
                end
                obj.solver_exists = true;
            end

            stats = obj.discrete_time_problem.solve();
            obj.stats = [obj.stats,stats];
            obj.w_all = [obj.w_all,obj.discrete_time_problem.w.res];
        end

        function clear_history(obj)
            obj.w_all = []; % Clear simulation data.
            obj.stats = []; % Clear simulation stats.
        end

        function [t_grid,x_res,t_grid_full,x_res_full] = simulate(obj, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            opts = obj.opts;
            x_res = obj.model.x0;
            x_res_full = obj.model.x0;
            t_grid = 0;
            t_grid_full = 0;
            obj.set_x0(obj.model.x0);
            obj.clear_history();
            t_current = 0;
            w0 = obj.discrete_time_problem.w.init;
            rbp = ~opts.right_boundary_point_explicit;

            for ii=1:opts.N_sim
                solver_stats = obj.solve(plugin);
                t_current = t_current + opts.T;
                if solver_stats.converged == 0
                    disp(['integrator_fesd: did not converge in step ', num2str(ii), ' constraint violation: ', num2str(solver_stats.constraint_violation, '%.2e')])
                elseif opts.print_level >=2
                    fprintf('Integration step %d / %d (%2.3f s / %2.3f s) converged in %2.3f s. \n',...
                        ii, opts.N_sim, t_current, opts.T_sim, solver_stats.cpu_time_total);
                end
                
                if obj.opts.dcs_mode ~= DcsMode.CLS
                    if rbp
                        x_step = obj.discrete_time_problem.w.x(0,0,opts.n_s).res;
                    else
                        x_step = [];
                    end
                    x_step = [x_step, obj.discrete_time_problem.w.x(:,:,opts.n_s+rbp).res];
                    x_step_full = obj.discrete_time_problem.w.x(:,:,:).res;
                    x_res = [x_res, x_step(:,2:end)];
                    x_res_full = [x_res_full, x_step_full(:,2:end)];
                    if opts.use_fesd
                        h = obj.discrete_time_problem.w.h(:,:).res;
                    else
                        h = ones(1,opts.N_finite_elements) * obj.discrete_time_problem.p.T().val/opts.N_finite_elements;
                    end
                    t_grid = [t_grid, t_grid(end) + cumsum(h)];
                    for ii = 1:length(h)
                        start = t_grid_full(end);
                        for jj = 1:opts.n_s
                            t_grid_full = [t_grid_full, start + opts.c_rk(jj)*h(ii)];
                        end
                        if rbp
                            t_grid_full = [t_grid_full, start + h(ii)];
                        end
                    end
                else
                    if rbp
                        x_step = obj.discrete_time_problem.w.x(0,0,opts.n_s).res;
                    else
                        x_step = [];
                    end
                    x_step = [x_step, obj.discrete_time_problem.w.x(:,:,[0,opts.n_s+rbp]).res];
                    x_step_full = obj.discrete_time_problem.w.x(:,:,:).res;
                    x_res = [x_res, x_step(:,(2+opts.no_initial_impacts):end)];
                    x_res_full = [x_res_full, x_step_full(:,(2+opts.no_initial_impacts):end)];
                    h = obj.discrete_time_problem.w.h(:,:).res;
                    for ii=1:length(h)
                        hi = h(ii);
                        if opts.no_initial_impacts && ii==1 
                            t_grid = [t_grid, t_grid(end)+hi];
                        else
                            t_grid = [t_grid, t_grid(end), t_grid(end)+hi];
                        end
                    end
                    for ii = 1:length(h)
                        start = t_grid_full(end);
                        for jj = 1:opts.n_s
                            t_grid_full = [t_grid_full, start + opts.c_rk(jj)*h(ii)];
                        end
                        if rbp
                            t_grid_full = [t_grid_full, start + h(ii)];
                        end
                    end
                end
                if opts.use_previous_solution_as_initial_guess
                    obj.discrete_time_problem.w.init = obj.discrete_time_problem.w.res;
                end
                obj.set_x0(x_step(:,end));
            end
        end

        function ret = get(obj, field)
            if isempty(obj.w_all)
                ret = [];
                return
            end

            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this integrator.']);
            end

            obj.discrete_time_problem.w.res = obj.w_all(:,1);
            if var.depth == 3
                rbp = ~opts.right_boundary_point_explicit;
                if opts.right_boundary_point_explicit
                    ret = var(:,:,obj.opts.n_s).res;
                else
                    ret = [var(0,0,obj.opts.n_s).res,...
                        var(1:opts.N_stages,1:opts.N_finite_elements(1),end).res];
                end
            else
                indexing(1:var.depth) = {':'};
                ret = var(indexing{:}).res;
            end

            for ii=2:size(obj.w_all, 2)
                obj.discrete_time_problem.w.res = obj.w_all(:,ii);
                if var.depth == 3
                    if opts.right_boundary_point_explicit
                        ret = [ret, var(1,:,obj.opts.n_s).res];
                    else
                        ret = [ret, var(1,1:opts.N_finite_elements(1),end).res];
                    end
                else
                    indexing(2:var.depth) = {':'};
                    if var.depth
                        indexing(1) = {1};
                    end
                    ret = [ret,var(indexing{:}).res];
                end 
            end
        end

        function ret = get_full(obj, field)
            if isempty(obj.w_all)
                ret = [];
                return
            end

            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this integrator.']);
            end

            obj.discrete_time_problem.w.res = obj.w_all(:,1);
            indexing(2:var.depth) = {':'};
            if var.depth
                indexing(1) = {':'};
            end
            ret = var(indexing{:}).res;

            for ii=2:size(obj.w_all, 2)
                obj.discrete_time_problem.w.res = obj.w_all(:,ii);
                indexing(2:var.depth) = {':'};
                if var.depth
                    indexing(1) = {1};
                end
                ret = [ret,var(indexing{:}).res];
            end
        end

        function t_grid = get_time_grid(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.get('h');
            else
                h = [];
                for ii=1:size(obj.w_all,2)
                    h = [h,ones(1,opts.N_finite_elements) * obj.discrete_time_problem.p.T().val/opts.N_finite_elements];
                end
            end
            t_grid = cumsum([0, h]);
        end

        function t_grid_full = get_time_grid_full(obj)
            opts = obj.opts;
            rbp = ~opts.right_boundary_point_explicit;
            if opts.use_fesd
                h = obj.get('h');
            else
                h = [];
                for ii=1:size(obj.w_all,2)
                    h = [h,ones(1,opts.N_finite_elements) * obj.discrete_time_problem.p.T().val/opts.N_finite_elements];
                end
            end
            t_grid = cumsum([0, h]);
            t_grid_full = 0;
            for ii = 1:length(h)
                start = t_grid_full(end);
                for jj = 1:opts.n_s
                    t_grid_full = [t_grid_full; start + opts.c_rk(jj)*h(ii)];
                end
                if rbp
                    t_grid_full = [t_grid_full; start + h(ii)];
                end
            end
        end

        function set(obj, varname, field, indices, value)
            if ~obj.discrete_time_problem.w.has_var(varname)
                error(['nosnoc:' char(varname) ' is not a valid field for this integrator.']);
            end
            var = obj.discrete_time_problem.w.(varname);
            var(indices{:}).(field) = value;
        end

        function set_x0(obj, x0)
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).init = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).lb = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).ub = x0;
        end

        function set_param(obj, param, value)
        % TODO (@anton) figure out how to do a set with indexing
            if ~obj.discrete_time_problem.p.has_var(param);
                error(['nosnoc:' char(param) ' does not exist as a parameter for this OCP']);
            end
            obj.discrete_time_problem.p.(param)().val = value;
        end

    end
end
