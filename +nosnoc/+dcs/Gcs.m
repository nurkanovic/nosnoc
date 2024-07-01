classdef Gcs < nosnoc.dcs.Base
    properties
        lambda %casadi.SX or casadi.MX: Lagrange multipliers for the gradient complementarity system.
        c_lift %casadi.SX or casadi.MX: Lift variables for the gap functions in :class:`nosnoc.model.Pds`.

        f_x % casadi.SX or casadi.MX: Right hand side of the gradient complementarity system

        dims % struct: dimensions struct.

        f_unconstrained_fun % casadi.Function: Unconstrained dynamics
        nabla_c_fun % casadi.Function: Gradient of the gap functions.
        c_fun % casadi.Function: Gap Functions
        g_c_lift_fun % casadi.Function: Lifting constraints for the gap functions in :class:`nosnoc.model.Pds`.
    end

    methods
        function obj = Gcs(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            dims = obj.dims;
            % dimensions
            dims.n_lambda = dims.n_c;
            if opts.gcs_lift_gap_functions
                dims.n_c_lift = dims.n_c;
            else
                dims.n_c_lift = 0;
            end
            obj.lambda = define_casadi_symbolic(opts.casadi_symbolic_mode,'lambda',dims.n_lambda);
            obj.c_lift = define_casadi_symbolic(opts.casadi_symbolic_mode,'c_lift',dims.n_c_lift);
            
            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            nabla_c = model.c.jacobian(model.x)';
            if opts.gcs_lift_gap_functions
                g_c_lift = obj.c_lift - model.c;
                c = obj.c_lift;
            else
                g_c_lift = [];
                c = model.c;
            end
            
            obj.f_x = model.f_x_unconstrained + model.E*nabla_c*obj.lambda;

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {obj.f_x});
            obj.f_unconstrained_fun = Function('f_unconstrained', {model.x, model.z, model.u, model.v_global, model.p}, {model.f_x_unconstrained});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {[]});
            obj.c_fun = Function('c_fun', {model.x, obj.c_lift, model.z, model.v_global, model.p}, {c});
            obj.g_c_lift_fun = Function('c_fun', {model.x, obj.c_lift, model.z, model.v_global, model.p}, {g_c_lift});
            obj.nabla_c_fun = Function('c_fun', {model.x, model.z, model.v_global, model.p}, {nabla_c});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.G_path_fun  = Function('G_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.G_path});
            obj.H_path_fun  = Function('H_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.H_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});
        end
    end

    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            propgrp = getPropertyGroups@nosnoc.dcs.Base(obj);
            group_title = 'Variables';
            var_list = struct;
            var_list.x = obj.model.x;
            if ~any(size(obj.model.u) == 0)
                var_list.u = obj.model.u;
            end
            if ~any(size(obj.model.z) == 0)
                var_list.z = obj.model.z;
            end
            var_list.lambda = obj.lambda;
            if ~any(size(obj.c_lift) == 0)
                var_list.c_lift = obj.c_lift;
            end
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end
end
