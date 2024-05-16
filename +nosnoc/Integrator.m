classdef Integrator < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_time_problem
        stats
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
                N_fe = sum(opts.N_finite_elements)
                opts.N_finite_elements = N_fe;
                opts.N_stages = 1;
            end
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
                elseif opts.dcs_mode == DcsMode.Step % TODO: RENAME
                    error("not implemented")
                else
                    error("PSS models can only be reformulated using the Stewart or Heaviside Step reformulations.")
                end
              case "nosnoc.model.heaviside"
                error("not implemented")
              case "nosnoc.model.cls"
                error("not implemented")
              case "nosnoc.model.pds"
                error("not implemented")
            end
        end

        function solve(obj, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);

            obj.stats = obj.discrete_time_problem.solve();
        end

        function [t_grid,x_res] = simulate(obj)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            opts = obj.opts;
            x_res = obj.model.x0;
            t_grid = 0;
            obj.set_x0(obj.model.x0);

            for ii=1:opts.N_sim
                obj.solve(plugin);
                x_step = obj.get_x();
                x_res = [x_res, x_step(:,2:end)];
                h = obj.discrete_time_problem.w.h(:,:).res;
                t_grid = [t_grid, t_grid(end) + cumsum(h)];
                obj.set_x0(x_step(:,end));
            end
        end

        function t_grid = get_time_grid(obj)
            h = obj.discrete_time_problem.w.h(:,:).res;
            t_grid = cumsum([0, h]);
        end

        function x = get_x(obj)
            opts = obj.opts;
            if opts.right_boundary_point_explicit
                x = obj.discrete_time_problem.w.x(:,:,obj.opts.n_s).res;
            else
                x = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
            end
        end

        function x = get_xend(obj)
            if opts.right_boundary_point_explicit
                x = obj.discrete_time_problem.w.x(1,opts.N_finite_elements,obj.opts.n_s).res;
            else
                x = obj.discrete_time_problem.w.x(1,opts.N_finite_elements,obj.opts.n_s+1).res;
            end
        end

        function set_x0(obj, x0)
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).init = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).lb = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).ub = x0;
        end
    end
end