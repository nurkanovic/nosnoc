classdef Gcs < nosnoc.dcs.Base
    properties
        lambda

        f_x

        dims

        nabla_c_fun
        c_fun
        g_comp_path_fun
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
            obj.lambda = define_casadi_symbolic(opts.casadi_symbolic_mode,'lambda',dims.n_lambda);

            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            nabla_c = model.c.jacobian(model.x)';
            
            obj.f_x = model.f_x + model.E*nabla_c*obj.lambda;

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.lambda, model.u, model.v_global, model.p}, {[]});
            obj.c_fun = Function('c_fun', {model.x, model.z, model.v_global, model.p}, {model.c});
            obj.nabla_c_fun = Function('c_fun', {model.x, model.z, model.v_global, model.p}, {nabla_c});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.g_comp_path_fun  = Function('g_comp_path', {model.x, model.z, model.u, model.v_global, model.p}, {[model.G_path,model.H_path]});
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
            var_list.lambda;
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end
end
