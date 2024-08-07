classdef Solver < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_time_problem
        stats
    end

    methods
        function obj = Solver(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            % Always process model and options
            opts.preprocess();
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
                    error("nosnoc: PSS models can only be reformulated using the Stewart or Heaviside Step reformulations.")
                end
              case "nosnoc.model.Heaviside"
                obj.dcs = nosnoc.dcs.Heaviside(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Heaviside(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Cls"
                if ~opts.use_fesd
                    error("nosnoc: The FESD-J reformulation only makes sense with use_fesd=true.")
                end
                obj.dcs = nosnoc.dcs.Cls(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Cls(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Pds"
                if ~opts.right_boundary_point_explicit
                    error("nosnoc: You are using an rk scheme with its right boundary point (c_n) not equal to one. Please choose another scheme e.g. RADAU_IIA.")
                end
                if opts.rk_representation == RKRepresentation.differential
                    error("nosnoc: Differential representation without lifting is unsupported for gradient complementarity systems. Use integral or lifted differential representation.")
                end
                obj.dcs = nosnoc.dcs.Gcs(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Gcs(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              otherwise
                error("nosnoc: Unknown model type.")
            end
        end

        function solve(obj, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);

            obj.stats = obj.discrete_time_problem.solve();
        end

        function ret = get(obj, field)
            opts = obj.opts;
            try
                warning off vdx:indexing:dot_reference_returns_vdx_var
                var = obj.discrete_time_problem.w.(field);
                warning on vdx:indexing:dot_reference_returns_vdx_var
            catch
                if strcmp(field, 'T_final')
                    warning("You are trying to get the final time from a non-time-optimal problem. Instead returning the fixed time.")
                    ret = obj.discrete_time_problem.p.T.val;
                    return
                else
                    error(['nosnoc:' char(field) ' is not a valid field for this OCP.']);
                end
                % TODO@anton print list of valid fields.
            end
            if var.depth == 3
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
        end

        function set_param(obj, param, value)
        % TODO (@anton) figure out how to do a set with indexing
            if ~obj.discrete_time_problem.p.has_var(param);
                error(['nosnoc:' char(param) ' does not exist as a parameter for this OCP.']);
            end
            warning off vdx:indexing:dot_reference_returns_vdx_var
            obj.discrete_time_problem.p.(param)().val = value;
            warning on vdx:indexing:dot_reference_returns_vdx_var
        end

        function ret = get_full(obj, field)
            opts = obj.opts;
            try
                warning off vdx:indexing:dot_reference_returns_vdx_var
                var = obj.discrete_time_problem.w.(field);
                warning on vdx:indexing:dot_reference_returns_vdx_var
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this OCP.']);
                % TODO @anton print list of valid fields.
            end
            indexing(1:var.depth) = {':'};
            ret = var(indexing{:}).res;
        end

        function t_grid = get_time_grid(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(opts.N_finite_elements))*(ones(1, sum(opts.N_finite_elements)));
            end
            if opts.use_speed_of_time_variables
                if opts.local_speed_of_time_variable
                    sot = repelem(obj.get("sot"), opts.N_finite_elements);
                    h = sot.*h;
                else
                    sot = obj.get("sot");
                    h = sot*h;
                end
            end
            t_grid = cumsum([0, h]);
        end

        function t_grid_full = get_time_grid_full(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(opts.N_finite_elements))*(ones(1, sum(opts.N_finite_elements)));
            end
            t_grid_full = 0;
            if opts.use_speed_of_time_variables
                if opts.local_speed_of_time_variable
                    sot = repelem(obj.get("sot"), opts.N_finite_elements);
                    h = sot.*h;
                else
                    sot = obj.get("sot");
                    h = sot*h;
                end
            end
            for ii = 1:length(h)
                start = t_grid_full(end);
                for jj = 1:opts.n_s
                    t_grid_full = [t_grid_full; start + opts.c_rk(jj)*h(ii)];
                end
            end
        end

        function t_grid = get_control_grid(obj)
            t_grid = [0];
            for ii=1:obj.opts.N_stages
                if obj.opts.use_fesd
                    h_sum = sum(obj.discrete_time_problem.w.h(ii,:).res);
                else
                    h_sum = obj.discrete_time_problem.p.T().val/obj.opts.N_stages;
                end
                t_grid = [t_grid, t_grid(end)+h_sum];
            end
        end

        function set(obj, varname, field, indices, value)
            if ~obj.discrete_time_problem.w.has_var(varname)
                error(['nosnoc:' char(varname) ' is not a valid field for this integrator.']);
            end
            var = obj.discrete_time_problem.w.(varname);
            var(indices{:}).(field) = value;
        end
    end
end
