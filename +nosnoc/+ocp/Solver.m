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
                error("not implemented")
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
                if strcmp(field, 'T_final')
                    warning("You are trying to get the final time from a non-time-optimal problem. Instead returning the fixed time.")
                    ret = obj.discrete_time_problem.p.T.val;
                    return
                else
                    error(['nosnoc:' char(field) ' is not a valid field for this OCP']);
                end
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

        function ret = get_full(obj, field)
            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error(['nosnoc:' char(field) ' is not a valid field for this OCP']);
                % TODO@anton print list of valid fields.
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
