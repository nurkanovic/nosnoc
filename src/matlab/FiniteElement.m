classdef FiniteElement < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_v
        ind_theta
        ind_lam
        ind_mu
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_h
        ind_elastic
        ind_boundary % index of bundary value lambda and mu, TODO is this even necessary?

        prev_fe
    end

    methods
        function obj = FiniteElement(prev_fe, settings, model, ctrl_idx, fe_idx)
            obj@NosnocFormulationObject();
            obj.ind_x = cell(settings.n_s, 1);
            obj.ind_v = cell(settings.n_s, 1);
            obj.ind_theta = cell(settings.n_s, settings.n_simplex);
            obj.ind_lam = cell(settings.n_s, settings.n_simplex);
            obj.ind_mu = cell(settings.n_s, settings.n_simplex);
            obj.ind_alpha = cell(settings.n_s, settings.n_simplex);
            obj.ind_lambda_n = cell(settings.n_s, settings.n_simplex);
            obj.ind_lambda_p = cell(settings.n_s, settings.n_simplex);
            obj.ind_h = [];
            obj.ind_elastic = [];
            obj.ind_boundary = [];

            obj.prev_fe = prev_fe;

            % TODO: create the finite element
            h = SX.sym(['h_' num2str(ctrl_idx) '_' num2str(fe_idx)]);
            h0 = h_ctrl_stage / model.N_finite_elements;
            ubh = (1 + settings.gamma_h) * h0;
            lbh = (1 - settings.gamma_h) * h0;
            obj.addVariable(h, 'h', lbh, ubh, h0);

            % RK stage stuff
            for ii = 1:settings.n_s
                % state / state derivative variables
                if settings.irk_representation == IrkRepresentation.DIFFERENTIAL
                    obj.addVariable(SX.sym(['V_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)'], model.nx),...
                                    'v',...
                                    -inf * ones(model.nx),...
                                    inf * ones(model.nx),...
                                    zeros(model.nx),...
                                    ii);
                end
                if settings.irk_representation == IrkRepresentation.INTEGRAL || settings.lift_irk_differential
                    obj.addVariable(SX.sym(['X_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)], model.nx),...
                                    'x',...
                                    -inf * ones(model.nx),...
                                    inf * ones(model.nx),...
                                    model.x0,...
                                    ii);
                end
                % algebraic variables
                if settings.pss_mode == PssMode.STEWART
                    % add thetas
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['theta_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], model.n_f_sys(ij)),...
                                        'theta',...
                                        zeros(model.n_f_sys(ij)),...
                                        inf * ones(model.n_f_sys(ij)),...
                                        settings.init_theta * ones(model.n_f_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambdas
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], model.n_f_sys(ij)),...
                                        'lam',...
                                        zeros(model.n_f_sys(ij)),...
                                        inf * ones(model.n_f_sys(ij)),...
                                        settings.init_lambda * ones(model.n_f_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add mu
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['mu_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], 1),...
                                        'mu',...
                                        -inf * ones(1),...
                                        inf * ones(1),...
                                        settings.init_mu * ones(1),...
                                        ii,...
                                        ij);
                    end
                elseif settings.pss_mode == PssMode.STEP
                    % add alpha
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['alpha_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], model.n_c_sys(ij)),...
                                        'alpha',...
                                        zeros(model.n_c_sys(ij)),...
                                        ones(model.n_c_sys(ij)),...
                                        settings.init_theta * ones(model.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_n
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], model.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(model.n_c_sys(ij)),...
                                        inf * ones(model.n_c_sys(ij)),...
                                        settings.init_lambda * ones(model.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)],model.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(model.n_c_sys(ij)),...
                                        inf * ones(model.n_c_sys(ij)),...
                                        settings.init_mu * ones(model.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                end
            end
            % Add right boundary points if needed
            if create_right_boundary_point
                if settings.pss_mode == PssMode.STEWART
                    % add lambdas
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], model.n_f_sys(ij)),...
                                        'lam',...
                                        zeros(model.n_f_sys(ij)),...
                                        inf * ones(model.n_f_sys(ij)),...
                                        settings.init_lambda * ones(model.n_f_sys(ij)),...
                                        settings.n_s,...
                                        ij);
                    end
                    % add mu
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['mu_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], 1),...
                                        'mu',...
                                        -inf * ones(1),...
                                        inf * ones(1),...
                                        settings.init_mu * ones(1),...
                                        settings.n_s,...
                                        ij);
                    end
                elseif settings.pss_mode == PssMode.STEP
                    % add lambda_n
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], model.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(model.n_c_sys(ij)),...
                                        inf * ones(model.n_c_sys(ij)),...
                                        settings.init_lambda * ones(model.n_c_sys(ij)),...
                                        settings.n_s,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], model.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(model.n_c_sys(ij)),...
                                        inf * ones(model.n_c_sys(ij)),...
                                        settings.init_mu * ones(model.n_c_sys(ij)),...
                                        settings.n_s,...
                                        ij);
                    end
                end
            end
            % add final X variables
            obj.addVariable(SX.sym(['X_end_' num2str(ctrl_idx) '_' num2str(fe_idx)], model.nx),...
                            'x',...
                            -inf * ones(model.nx),...
                            inf * ones(model.nx),...
                            model.x0,...
                            -1);
        end
        
        function lambda = lambda(obj, varargin)
            import casadi.*
            p = inputParser();
            p.FunctionName('lambda')
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'stage',[]);
            addOptional(p, 'sys',[]);
            parse(p, obj, symbolic, lb, ub, initial, idx, varargin{:});

            
            if ismember('stage', p.UsingDefaults)
                lambda = vertcat(obj.w([obj.ind_lam{:}]), obj.w([obj.ind_lambda_n{:}]), obj.w([obj.ind_lambda_p{:}]));
            else
                if ~ismember('sys', p.UsingDefaults)
                    lambda = vertcat(obj.w([obj.ind_lam{p.Results.stage, p.Results.sys}]),...
                                     obj.w([obj.ind_lambda_n{p.Results.stage, p.Results.sys}]),...
                                     obj.w([obj.ind_lambda_p{p.Results.stage, p.Results.sys}]));
                else
                    lambda = vertcat(obj.w([obj.ind_lam{p.Results.stage, :}]),...
                                     obj.w([obj.ind_lambda_n{p.Results.stage, :}]),...
                                     obj.w([obj.ind_lambda_p{p.Results.stage, :}]));
                end
            end
        end

        function sum_lambda = sumLambda(obj, varargin)
            import casadi.*
            p = inputParser();
            p.FunctionName('sumLambda')
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, symbolic, lb, ub, initial, idx, varargin{:});

            n_stages = shape(obj.ind_lam, 1);
            if ismember('sys', p.UsingDefaults)
                lambdas = obj.lambda(1:n_stages);
            else
                lambdas = obj.lambda(1:n_stages, p.Results.sys);
            end
            lambdas = [lambdas, obj.prev_fe.lambda()]
            sum_lambda = sum(lambdas, 2);
        end

        function z = rkStageZ(obj, stage)
            import casadi.*

            idx = [[obj.ind_theta{stage, :}],...
                   [ebj.ind_lam{stage, :}],...
                   [obj.ind_mu{stage, :}],...
                   [obj.ind_alpha{stage, :}],...
                   [obj.ind_lambda_n{stage, :}],...
                   [obj.ind_lambda_p{stage, :}]];

            z = obj.w(idx);
        end

        function theta = theta(obj, varargin)
            import casadi.*
            p = inputParser();
            p.FunctionName('theta')
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'stage',[]);
            addOptional(p, 'sys',[]);
            parse(p, obj, symbolic, lb, ub, initial, idx, varargin{:});
            
            if ismember('stage', p.UsingDefaults)
                n_alpha = length([obj.ind_alpha{:}]);
                theta = vertcat(obj.w([obj.ind_theta{:}]),...
                                obj.w([obj.ind_alpha{:}]),...
                                ones(n_alpha, 1) - obj.w([obj.ind_alpha{:}]));
            else
                if ~ismember('sys', p.UsingDefaults)
                    n_alpha = length([obj.ind_alpha{p.Results.stage, p.Results.sys}]);
                    theta = vertcat(obj.w([obj.ind_theta{p.Results.stage, p.Results.sys}]),...
                                    obj.w([obj.ind_alpha{p.Results.stage, p.Results.sys}]),...
                                    ones(n_alpha, 1) - obj.w([obj.ind_alpha{p.Results.stage, p.Results.sys}]));
                else
                    n_alpha = length([obj.ind_alpha{p.Results.stage, :}]);
                    theta = vertcat(obj.w([obj.ind_theta{p.Results.stage, :}]),...
                                    obj.w([obj.ind_alpha{p.Results.stage, :}]),...
                                    ones(n_alpha, 1) - obj.w([obj.ind_lambda_p{p.Results.stage, :}]));
                end
            end
        end

        function sum_theta = sumTheta(obj)

        end

        function forwardSimulation(obj)

        end

        function createComplementarityConstraints(obj)

        end

        function stepEquilibration(obj)

        end
    end
end

