classdef FiniteElementZero < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_lam
        ind_lambda_n
        ind_lambda_p

        ctrl_idx
        fe_idx

        prev_fe
    end

    properties(Dependent, SetAccess=private, Hidden)
        x
        v
        theta
        lam
        mu
        alpha
        lambda_n
        lambda_p
        h

        lambda
    end
    
    methods
        function obj = FiniteElementZero(settings, dims, model)
            import casadi.*
            obj@NosnocFormulationObject();

            obj.ind_x = cell(1, 1);
            obj.ind_lam = cell(1,dims.n_sys);
            obj.ind_lambda_n = cell(1,dims.n_sys);
            obj.ind_lambda_p = cell(1,dims.n_sys);

            % X0
            if settings.there_exist_free_x0
                x0_lb = model.x0;
                x0_ub = model.x0;
                x0_ub(model.ind_free_x0) = inf;
                x0_lb(model.ind_free_x0) = -inf;

                x0_ub(model.ind_free_x0) = ubx(model.ind_free_x0);
                x0_lb(model.ind_free_x0) = lbx(model.ind_free_x0);

                x0 = define_casadi_symbolic(settings.casadi_symbolic_mode, 'X0', dims.n_x);
                obj.addVariable(x0,...
                                'x',...
                                x0_lb,...
                                x0_ub,...
                                model.x0,...
                                1);
            else
                x0 = define_casadi_symbolic(settings.casadi_symbolic_mode, 'X0', dims.n_x);
                obj.addVariable(x0,...
                                'x',...
                                model.x0,...
                                model.x0,...
                                model.x0,...
                                1);
            end

            % lambda00
            % TODO Use define_casadi_symbolic
            if settings.pss_mode == PssMode.Stewart
                for ij=1:dims.n_sys
                    lam0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_' num2str(ij)], dims.n_f_sys(ij));
                    obj.addVariable(lam0,...
                                    'lam',...
                                    -inf * ones(dims.n_f_sys(ij),1),...
                                    inf * ones(dims.n_f_sys(ij),1),...
                                    settings.initial_lambda * ones(dims.n_f_sys(ij),1),...
                                    1,...
                                    ij);
                end
            elseif settings.pss_mode == PssMode.Step
                for ij=1:dims.n_sys
                    lambda_n0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_n_' num2str(ij)], dims.n_c_sys(ij));
                    lambda_p0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_p_' num2str(ij)], dims.n_c_sys(ij));
                    obj.addVariable(lambda_n0,...
                                    'lambda_n',...
                                    zeros(dims.n_c_sys(ij),1),...
                                    inf * ones(dims.n_c_sys(ij),1),...
                                    settings.initial_lambda_0 * ones(dims.n_c_sys(ij),1),...
                                    1,...
                                    ij);
                    obj.addVariable(lambda_p0,...
                                    'lambda_p',...
                                    zeros(dims.n_c_sys(ij),1),...
                                    inf * ones(dims.n_c_sys(ij),1),...
                                    settings.initial_lambda_1 * ones(dims.n_c_sys(ij),1),...
                                    1,...
                                    ij);
                end
            end
        end

        function lambda = get.lambda(obj)
            import casadi.*
            grab = @(l, ln, lp) vertcat(obj.w(l), obj.w(ln), obj.w(lp));

            lambda = cellfun(grab, obj.ind_lam, obj.ind_lambda_p, obj.ind_lambda_n, 'UniformOutput', false);
        end

        function x = get.x(obj)
            x= cellfun(@(x) obj.w(x), obj.ind_x, 'UniformOutput', false);
        end
    end
end
