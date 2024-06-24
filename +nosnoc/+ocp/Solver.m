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
                elseif opts.dcs_mode == DcsMode.Heaviside % TODO: RENAME
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
                obj.dcs = nosnoc.dcs.Gcs(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Gcs(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              otherwise
                error("Unknown model type")
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
                var = obj.discrete_time_problem.w.(field);
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this OCP']);
                % TODO@anton print list of valid fields.
            end
            if var.depth == 3
                if opts.right_boundary_point_explicit
                    ret = var(:,:,obj.opts.n_s).res;
                else
                    ret = [var(0,0,obj.opts.n_s).res,...
                        var(1:opts.N_stages,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
                end
            else
                indexing(1:var.depth) = {':'};
                ret = var(indexing{:}).res;
            end
        end

        function set_param(obj, param, value)
        % TODO (@anton) figure out how to do a set with indexing
            if ~obj.discrete_time_problem.p.has_var(param);
                error(['nosnoc:' char(param) ' does not exist as a parameter for this OCP']);
            end
            obj.discrete_time_problem.p.(param)().val = value;
        end

        function ret = get_full(obj, field)
            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this OCP']);
                % TODO @anton print list of valid fields.
            end
            indexing(1:var.depth) = {':'};
            ret = var(indexing{:}).res;
        end

        function t_grid = get_time_grid(obj)
            if obj.opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(obj.opts.N_finite_elements))*(ones(1, sum(obj.opts.N_finite_elements)));
            end
            t_grid = cumsum([0, h]);
        end

        function t_grid_full = get_time_grid_full(obj)
            if obj.opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(obj.opts.N_finite_elements))*(ones(1, sum(obj.opts.N_finite_elements)));
            end
            t_grid_full = 0;
            for ii = 1:length(h)
                for jj = 1:opts.n_s
                    t_grid_full = [t_grid_full; t_grid_full(end) + opts.c_rk(jj)*h(ii)];
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
    end
end
