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
            h_ctrl_stage = settings.terminal_time / settings.N_stages;
            h0 = array([h_ctrl_stage / array(settings.Nfe_list[ctrl_idx])]);
            ubh = (1 + settings.gamma_h) * h0;
            lbh = (1 - settings.gamma_h) * h0;
            obj.addVariable(h, 'h', lbh, ubh, h0);

            % RK stage stuff
            for ii = 1:settings.n_s
                % state / state derivative variables
                if settings.irk_representation == IrkRepresentation.DIFFERENTIAL
                    obj.addVariable(SX.sym(['V_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)'], model.n_x),...
                                    'v',...
                                    -inf * ones(model.n_x),...
                                    inf * ones(model.n_x),...
                                    zeros(model.n_x),...
                                    ii);
                end
                if settings.irk_representation == IrkRepresentation.INTEGRAL || settings.lift_irk_differential
                    obj.addVariable(SX.sym(['X_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)], model.n_x),...
                                    'x',...
                                    -inf * ones(model.n_x),...
                                    inf * ones(model.n_x),...
                                    model.x0,...
                                    ii);
                end
                % algebraic variables
                if settings.pss_mode == PssMode.STEWART
                    % add thetas
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['theta_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'theta',...
                                        zeros(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.init_theta * ones(dims.n_f_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambdas
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'lam',...
                                        zeros(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.init_lambda * ones(dims.n_f_sys(ij)),...
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
                        obj.addVariable(SX.sym(['alpha_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'alpha',...
                                        zeros(dims.n_c_sys(ij)),...
                                        ones(dims.n_c_sys(ij)),...
                                        settings.init_theta * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_n
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.init_lambda * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)],dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.init_mu * ones(dims.n_c_sys(ij)),...
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
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'lam',...
                                        zeros(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.init_lambda * ones(dims.n_f_sys(ij)),...
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
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.init_lambda * ones(dims.n_c_sys(ij)),...
                                        settings.n_s,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:model.n_simplex
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.init_mu * ones(dims.n_c_sys(ij)),...
                                        settings.n_s,...
                                        ij);
                    end
                end
            end
            % add final X variables
            obj.addVariable(SX.sym(['X_end_' num2str(ctrl_idx) '_' num2str(fe_idx)], model.n_x),...
                            'x',...
                            -inf * ones(model.n_x),...
                            inf * ones(model.n_x),...
                            model.x0,...
                            -1);
        end
        
        function lambda = Lambda(obj, varargin)
            
        end
    end
end
