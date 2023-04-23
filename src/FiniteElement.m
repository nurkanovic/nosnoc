% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.
classdef FiniteElement < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_v
        ind_z
        % Stewart
        ind_theta
        ind_lam
        ind_mu
        % Step
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_beta
        ind_theta_step
        % CLS
        ind_x_left_bp
        ind_lambda_normal
        ind_lambda_tangent
        ind_y_gap
        % friction multipliers and lifting
        % conic
        ind_gamma
        ind_beta_conic
        % poly
        ind_gamma_d
        ind_beta_d
        ind_delta_d
        % variables related to conic
        ind_p_vt
        ind_n_vt
        ind_alpha_vt
        % variables only at element boundary
        ind_Y_gap
        ind_Lambda_normal
        ind_Lambda_tangent
        %
        ind_Gamma
        ind_Beta_conic
        %
        ind_Gamma_d
        ind_Beta_d
        ind_Delta_d
        %
        ind_P_vn
        ind_N_vn
        ind_P_vt
        ind_N_vt
        ind_Alpha_vt
        % misc
        ind_h
        ind_nu_lift
        ind_elastic
        ind_boundary % index of bundary value lambda and mu, TODO is this even necessary?


        ind_eq
        ind_ineq
        ind_comp

        cross_comp_pairs

        ctrl_idx
        fe_idx

        model
        settings
        dims

        u

        T_final

        prev_fe
    end

    properties(Dependent, SetAccess=private, Hidden)
        x
        v
        % Stewart
        lam
        mu
        % Step
        alpha
        lambda_n
        lambda_p
        % CLS
        lambda_normal
        lambda_tangent
        y_gap
        % conic
        gamma
        beta_conic
        % poly
        gamma_d
        beta_d
        delta_d
        % variables related to conic
        p_vt
        n_vt
        alpha_vt
        % variables only at element boundary
        Y_Gap
        Lambda_normal
        Lambda_tangent
        Gamma
        Gamma_d
        Beta_conic
        Delta_d
        P_vn
        N_vn
        P_vt
        N_vt
        Alpha_vt
        % misc
        nu_lift
        h

        elastic

        nu_vector
    end

    properties(SetAccess=private)
        % TODO: rename into cross_comp_cont cross_comp_discont?
        cross_comp_cont_0 % cross comp variables that are cont (requite conecetion to prev. element, e.g. lambda in stewart)
        cross_comp_discont_0 % cross comp variables that are discont (e.g. theta in stewart)
        cross_comp_cont_1 % cross comp variables that are cont and get complemented with cross_comp_cont_2 that are cont as well
        cross_comp_cont_2
    end

    methods
        function obj = FiniteElement(prev_fe, settings, model, dims, ctrl_idx, fe_idx, varargin)
            import casadi.*
            obj@NosnocFormulationObject();

            p = inputParser();
            p.FunctionName = 'FiniteElement';

            % TODO: add checks.
            addRequired(p, 'prev_fe');
            addRequired(p, 'settings');
            addRequired(p, 'model');
            addRequired(p, 'dims');
            addRequired(p, 'ctrl_idx');
            addRequired(p, 'fe_idx');
            addOptional(p, 'T_final',[]);
            parse(p, prev_fe, settings, model, dims, ctrl_idx, fe_idx, varargin{:});

            if settings.right_boundary_point_explicit || settings.dcs_mode == 'CLS'
                rbp_allowance = 0;
            else
                rbp_allowance = 1;
            end

            if settings.dcs_mode == "CLS"
                obj.ind_x = cell(dims.n_s, 1);
            else
                obj.ind_x = cell(dims.n_s+rbp_allowance, 1);
            end
            obj.ind_v = cell(dims.n_s, 1);
            obj.ind_z = cell(dims.n_s+rbp_allowance, 1);
            % Stewart
            obj.ind_theta = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lam = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_mu = cell(dims.n_s+rbp_allowance, dims.n_sys);
            % Step
            obj.ind_alpha = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_n = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_p = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_theta_step = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_beta = cell(dims.n_s+rbp_allowance, 1);
            % CLS
            obj.ind_lambda_normal = cell(dims.n_s,1);
            obj.ind_lambda_tangent = cell(dims.n_s,1);
            obj.ind_y_gap = cell(dims.n_s,1);
            % friction multipliers and lifting
            % conic
            obj.ind_gamma = cell(dims.n_s,1);
            obj.ind_beta_conic = cell(dims.n_s,1);
            % poly
            obj.ind_gamma_d = cell(dims.n_s,1);
            obj.ind_beta_d = cell(dims.n_s,1);
            obj.ind_delta_d = cell(dims.n_s,1);
            % variables related to conic
            obj.ind_p_vt = cell(dims.n_s,1);
            obj.ind_n_vt = cell(dims.n_s,1);
            obj.ind_alpha_vt = cell(dims.n_s,1);
            % variables only at element boundary
            obj.ind_Lambda_normal = cell(1,1);
            obj.ind_Lambda_tangent = cell(1,1);
            obj.ind_Gamma = cell(1,1);
            obj.ind_Beta_d = cell(1,1);
            obj.ind_Gamma_d = cell(1,1);
            obj.ind_Beta_conic = cell(1,1);
            obj.ind_Delta_d = cell(1,1);
            obj.ind_P_vn = cell(1,1);
            obj.ind_N_vn = cell(1,1);
            obj.ind_P_vt = cell(1,1);
            obj.ind_N_vt = cell(1,1);
            obj.ind_Alpha_vt = cell(1,1);

            % misc
            obj.ind_h = [];
            obj.ind_elastic = [];
            obj.ind_boundary = [];

            obj.ctrl_idx = ctrl_idx;
            obj.fe_idx = fe_idx;

            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;

            obj.prev_fe = prev_fe;

            if settings.use_fesd
                h = define_casadi_symbolic(settings.casadi_symbolic_mode, ['h_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)]);
                h_ctrl_stage = model.T/dims.N_stages;
                h0 = h_ctrl_stage / dims.N_finite_elements(ctrl_idx);
                ubh = (1 + settings.gamma_h) * h0;
                lbh = (1 - settings.gamma_h) * h0;
                if settings.time_rescaling && ~settings.use_speed_of_time_variables
                    % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                    ubh = (1+settings.gamma_h)*h0*settings.s_sot_max;
                    lbh = (1-settings.gamma_h)*h0/settings.s_sot_min;
                end
                obj.addVariable(h, 'h', lbh, ubh, h0);
            end
            obj.T_final = p.Results.T_final;

            % Left boundary point needed for dcs_mode = cls (corresponding to t^+ (post impacts))
            if settings.dcs_mode == "CLS"
                x = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['X_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_0'],...
                    dims.n_x);
                obj.addVariable(x,...
                    'x_left_bp',...
                    model.lbx,...
                    model.ubx,...
                    model.x0,...
                    1);
            end
            % RK stage stuff
            for ii = 1:dims.n_s
                % state / state derivative variables
                if (settings.irk_representation == IrkRepresentation.differential ||...
                        settings.irk_representation == IrkRepresentation.differential_lift_x)
                    v = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['V_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)'],...
                        dims.n_x);
                    obj.addVariable(v,...
                        'v',...
                        -inf * ones(dims.n_x,1),...
                        inf * ones(dims.n_x,1),...
                        zeros(dims.n_x,1),...
                        ii);
                end
                if (settings.irk_representation == IrkRepresentation.integral ||...
                        settings.irk_representation == IrkRepresentation.differential_lift_x)
                    if settings.x_box_at_stg
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif settings.x_box_at_fe && ii == dims.n_s && settings.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif fe_idx == dims.N_finite_elements(ctrl_idx) && ii == dims.n_s && settings.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    else
                        lbx = -inf * ones(dims.n_x, 1);
                        ubx = inf * ones(dims.n_x, 1);
                    end
                    x = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['X_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_x);
                    obj.addVariable(x,...
                        'x',...
                        lbx,...
                        ubx,...
                        model.x0,...
                        ii);
                end
                % algebraic variables
                if settings.dcs_mode == DcsMode.Stewart
                    % add thetas
                    for ij = 1:dims.n_sys
                        theta = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['theta_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(theta,...
                            'theta',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_theta * ones(dims.n_f_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add lambdas
                    for ij = 1:dims.n_sys
                        lam = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_lambda * ones(dims.n_f_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf,...
                            inf,...
                            settings.initial_mu,...
                            ii,...
                            ij);
                    end
                elseif settings.dcs_mode == DcsMode.Step
                    % add alpha
                    for ij = 1:dims.n_sys
                        alpha = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['alpha_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(alpha,...
                            'alpha',...
                            zeros(dims.n_c_sys(ij), 1),...
                            ones(dims.n_c_sys(ij), 1),...
                            settings.initial_theta * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add lambda_n and lambda_p
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        lambda_p = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_0 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_1 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % TODO: Clean this up (maybe as a function to reduce the indent level.
                    if settings.time_freezing_inelastic
                        %                         if ~settings.pss_lift_step_functions
                        theta_step = define_casadi_symbolic(settings.casadi_symbolic_mode, ['theta_step_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_theta_step);
                        beta = define_casadi_symbolic(settings.casadi_symbolic_mode, ['beta_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_beta);
                        beta_guess = full(model.g_lift_beta_fun(settings.initial_alpha*ones(dims.n_alpha,1)));
                        theta_step_guess = full(model.g_lift_theta_step_fun(settings.initial_alpha*ones(dims.n_alpha,1),beta_guess));
                        obj.addVariable(theta_step,...
                            'theta_step',...
                            -inf*ones(dims.n_theta_step,1),...
                            inf*ones(dims.n_theta_step,1),...
                            theta_step_guess,...
                            ii);
                        if settings.pss_lift_step_functions
                            obj.addVariable(beta,...
                                'beta',...
                                -inf*ones(dims.n_beta,1),...
                                inf*ones(dims.n_beta,1),...
                                beta_guess,...
                                ii);
                        end
                    end
                elseif settings.dcs_mode == DcsMode.CLS
                    lambda_normal = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['lambda_normal_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_contacts);
                    obj.addVariable(lambda_normal,...
                        'lambda_normal',...
                        zeros(dims.n_contacts,1),...
                        inf * ones(dims.n_contacts, 1),...
                        ones(dims.n_contacts, 1),...
                        ii,1);
                    y_gap = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['y_gap_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_contacts);
                    obj.addVariable(y_gap,...
                        'y_gap',...
                        zeros(dims.n_contacts,1),...
                        inf * ones(dims.n_contacts, 1),...
                        ones(dims.n_contacts, 1),...
                        ii,1);
                    if model.friction_exists
                        lambda_tangent = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_tangent' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                            dims.n_tangents);

                        if isequal(settings.friction_model,'Polyhedral')
                            obj.addVariable(lambda_tangent,...
                                'lambda_tangent',...
                                -inf *ones(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);

                            gamma_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                ['gamma_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);
                            obj.addVariable(gamma_d,...
                                'gamma_d',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            beta_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                ['beta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);
                            obj.addVariable(beta_d,...
                                'beta_d',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);


                            delta_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                ['delta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_tangents);

                            obj.addVariable(delta_d,...
                                'delta_d',...
                                zeros(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);
                        end
                        if isequal(settings.friction_model,'Conic')
                            obj.addVariable(lambda_tangent,...
                                'lambda_tangent',...
                                zeros(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);

                            gamma = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                ['gamma' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);

                            obj.addVariable(gamma,...
                                'gamma',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            beta = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                ['beta' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);

                            obj.addVariable(beta,...
                                'beta',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            switch settings.conic_model_switch_handling
                                case 'Plain'
                                    % no extra vars
                                case 'Abs'
                                    p_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                        ['p_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(p_vt,...
                                        'p_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    n_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                        ['n_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(n_vt,...
                                        'n_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                case 'Lp'
                                    p_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                        ['p_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(p_vt,...
                                        'p_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    n_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                        ['n_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(n_vt,...
                                        'n_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    alpha_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                        ['alpha_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(alpha_vt,...
                                        'alpha_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                            end
                        end
                    end
                end
                % add user variables
                z = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['z_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                    dims.n_z);
                obj.addVariable(z,...
                    'z',...
                    model.lbz,...
                    model.ubz,...
                    model.z0,...
                    ii);
            end

            if settings.dcs_mode == DcsMode.CLS
                %  IMPULSE VARIABLES
                Lambda_normal = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['Lambda_normal_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    dims.n_contacts);
                obj.addVariable(Lambda_normal,...
                    'Lambda_normal',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    ones(dims.n_contacts, 1),1);
                P_vn = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['P_vn' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                    dims.n_contacts);
                obj.addVariable(P_vn,...
                    'P_vn',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    ones(dims.n_contacts, 1),1);
                N_vn = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['N_vn' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                    dims.n_contacts);
                obj.addVariable(N_vn,...
                    'N_vn',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    ones(dims.n_contacts, 1),1);

                Y_gap = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['Y_gap_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                    dims.n_contacts);
                obj.addVariable(Y_gap,...
                    'Y_gap',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    ones(dims.n_contacts, 1),1);
                if model.friction_exists
                    Lambda_tangent = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['Lambda_tangent' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                        dims.n_tangents);

                    if settings.friction_model == FrictionModel.Polyhedral
                        obj.addVariable(Lambda_tangent,...
                            'Lambda_tangent',...
                            -inf *ones(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);

                        Gamma_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['Gamma_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);
                        obj.addVariable(Gamma_d,...
                            'Gamma_d',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Beta_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['Beta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);
                        obj.addVariable(Beta_d,...
                            'Beta_d',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Delta_d = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['Delta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_tangents);

                        obj.addVariable(Delta_d,...
                            'Delta_d',...
                            zeros(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);
                    end
                    if settings.friction_model == FrictionModel.Conic
                        obj.addVariable(Lambda_tangent,...
                            'Lambda_tangent',...
                            zeros(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);

                        Gamma = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['Gamma' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);

                        obj.addVariable(Gamma,...
                            'Gamma',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Beta_conic = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['Beta_conic' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);

                        obj.addVariable(Beta_conic,...
                            'Beta_conic',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        switch settings.conic_model_switch_handling
                            case 'Plain'
                                % no extra vars
                            case 'Abs'
                                P_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                    ['P_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(P_vt,...
                                    'P_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                N_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                    ['N_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(N_vt,...
                                    'N_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                            case 'Lp'
                                P_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                    ['P_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(P_vt,...
                                    'P_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                N_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                    ['N_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(N_vt,...
                                    'N_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                Alpha_vt = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                                    ['Alpha_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(Alpha_vt,...
                                    'Alpha_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                        end
                    end
                end
            end


            % Add right boundary points if needed
            if ~settings.right_boundary_point_explicit
                if settings.dcs_mode == DcsMode.Stewart
                    % add lambdas
                    for ij = 1:dims.n_sys
                        lam = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_lambda * ones(dims.n_f_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf * ones(1),...
                            inf * ones(1),...
                            settings.initial_mu * ones(1),...
                            dims.n_s+1,...
                            ij);
                    end

                elseif settings.dcs_mode == DcsMode.Step
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_0 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add lambda_p
                    for ij = 1:dims.n_sys
                        lambda_p = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_1 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end
                end
            end

            if (~settings.right_boundary_point_explicit ||...
                    settings.irk_representation == IrkRepresentation.differential)
                if settings.x_box_at_stg || settings.x_box_at_fe || fe_idx == dims.N_finite_elements(ctrl_idx)
                    lbx = model.lbx;
                    ubx = model.ubx;
                else
                    lbx = -inf * ones(dims.n_x);
                    ubx = inf * ones(dims.n_x);
                end

                % add final X variables
                x = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['X_end_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    dims.n_x);
                obj.addVariable(x,...
                    'x',...
                    lbx,...
                    ubx,...
                    model.x0,...
                    dims.n_s+rbp_allowance);
            end
            if strcmpi(settings.step_equilibration, 'direct_homotopy_lift')
                nu_lift = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['nu_lift_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    1);
                obj.addVariable(nu_lift,...
                    'nu_lift',...
                    1,...
                    -inf,...
                    inf);
            end

            % calculate lambda and theta


            if model.friction_exists && settings.friction_model == 'Conic' && settings.conic_model_switch_handling == 'Abs'

                grab = @(l, ln, lp ,yg, g, pvt, nvt, gd,dd) vertcat(obj.w(l), obj.w(ln), obj.w(lp), obj.w(yg), obj.w(g), obj.w(gd), obj.w(dd));
                obj.cross_comp_cont_0 = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p,...
                    obj.ind_y_gap, obj.ind_gamma, obj.ind_p_vt, obj.ind_n_vt, obj.ind_gamma_d, obj.ind_delta_d, 'UniformOutput', false);

                grab = @(t, a, ln, b, avt, bd, lt) vertcat(obj.w(t), obj.w(a), ones(size(a))' - obj.w(a), obj.w(ln), obj.w(b), obj.w(bd), obj.w(lt));
                obj.cross_comp_discont_0 = cellfun(grab, obj.ind_theta, obj.ind_alpha,...
                    obj.ind_lambda_normal, obj.ind_beta_conic, obj.ind_alpha_vt, obj.ind_beta_d, obj.ind_lambda_tangent, 'UniformOutput', false);


                grab = @(pvt) vertcat(obj.w(pvt));
                obj.cross_comp_cont_1 = cellfun(grab, obj.ind_p_vt, 'UniformOutput', false);
                grab = @(nvt) vertcat(obj.w(nvt));
                obj.cross_comp_cont_2 = cellfun(grab,  obj.ind_n_vt, 'UniformOutput', false);
            else
                % ??
                grab = @(l, ln, lp ,yg, g, pvt, nvt, gd,dd) vertcat(obj.w(l), obj.w(ln), obj.w(lp), obj.w(yg), obj.w(g), obj.w(pvt), obj.w(nvt), obj.w(gd), obj.w(dd));
                obj.cross_comp_cont_0 = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p,...
                    obj.ind_y_gap, obj.ind_gamma, obj.ind_p_vt, obj.ind_n_vt, obj.ind_gamma_d, obj.ind_delta_d, 'UniformOutput', false);

                grab = @(t, a, ln, b, avt, bd, lt) vertcat(obj.w(t), obj.w(a), ones(size(a))' - obj.w(a), obj.w(ln), obj.w(b), ones(size(avt))' - obj.w(avt), obj.w(avt), obj.w(bd), obj.w(lt));
                obj.cross_comp_discont_0 = cellfun(grab, obj.ind_theta, obj.ind_alpha,...
                    obj.ind_lambda_normal, obj.ind_beta_conic, obj.ind_alpha_vt, obj.ind_beta_d, obj.ind_lambda_tangent, 'UniformOutput', false);
            end


        end

        function h = get.h(obj)
            if obj.settings.use_fesd
                h = obj.w(obj.ind_h);
            elseif obj.settings.time_optimal_problem && ~obj.settings.use_speed_of_time_variables
                h = obj.T_final/(obj.dims.N_stages*obj.dims.N_finite_elements(obj.ctrl_idx));
            else
                h = obj.model.T/(obj.dims.N_stages*obj.dims.N_finite_elements(obj.ctrl_idx));
            end
        end

        function x = get.x(obj)
            x = cellfun(@(x) obj.w(x), obj.ind_x, 'UniformOutput', false);
        end

        function v = get.v(obj)
            v = cellfun(@(v) obj.w(v), obj.ind_v, 'UniformOutput', false);
        end

        function elastic = get.elastic(obj)
            elastic = obj.w(obj.ind_elastic);
        end

        function nu_lift = get.nu_lift(obj)
            nu_lift = obj.w(obj.ind_nu_lift);
        end

        function nu_vector = get.nu_vector(obj)
            import casadi.*
            if obj.settings.use_fesd && obj.fe_idx > 1
                eta_k = obj.prev_fe.sumCrossCompCont0().*obj.sumCrossCompCont0() + obj.prev_fe.sumCrossCompDiscont0().*obj.sumCrossCompDiscont0();
                nu_vector = 1;
                for jjj=1:length(eta_k)
                    nu_vector = nu_vector * eta_k(jjj);
                end
            else
                nu_vector = [];
            end
        end

        function sum_cross_comp_cont_0 = sumCrossCompCont0(obj, varargin)
            import casadi.*
            p = inputParser();
            p.FunctionName = 'sumCrossCompCont0';

            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});

            if ismember('sys', p.UsingDefaults)
                cross_comp_cont_0 = obj.cross_comp_cont_0;
                cross_comp_cont_0_vec = arrayfun(@(row) vertcat(cross_comp_cont_0{row, :}), 1:size(cross_comp_cont_0,1), 'UniformOutput', false);
                cross_comp_cont_0_vec = [cross_comp_cont_0_vec, {vertcat(obj.prev_fe.cross_comp_cont_0{end,:})}];
            else
                cross_comp_cont_0_vec = obj.cross_comp_cont_0(:,p.Results.sys).';
                cross_comp_cont_0_vec = [cross_comp_cont_0_vec, obj.prev_fe.cross_comp_cont_0(end,p.Results.sys)];
            end
            sum_cross_comp_cont_0 = sum([cross_comp_cont_0_vec{:}], 2);
        end



        function sum_cross_comp_discont_0 = sumCrossCompDiscont0(obj, varargin)
            p = inputParser();
            p.FunctionName = 'sumCrossCompDiscont0';

            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});
            %obj.theta
            if ismember('sys', p.UsingDefaults)
                cross_comp_discont_0 = obj.cross_comp_discont_0;
                cross_comp_discont_0_vec = arrayfun(@(row) vertcat(cross_comp_discont_0{row, :}), 1:size(cross_comp_discont_0,1), 'UniformOutput', false);
            else
                cross_comp_discont_0_vec = obj.cross_comp_discont_0(:,p.Results.sys).';
            end
            sum_cross_comp_discont_0 = sum([cross_comp_discont_0_vec{:}], 2);
        end

        function sum_cross_comp_cont_1 = sumCrossCompCont1(obj, varargin)
            % TODO!!! : DO appropiate cross comp functions
            import casadi.*
            p = inputParser();
            p.FunctionName = 'sumCrossCompCont1';

            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});

            if ismember('sys', p.UsingDefaults)
                cross_comp_cont_1 = obj.cross_comp_cont_1;
                cross_comp_cont_1_vec = arrayfun(@(row) vertcat(cross_comp_cont_1{row, :}), 1:size(cross_comp_cont_1,1), 'UniformOutput', false);
                cross_comp_cont_1_vec = [cross_comp_cont_1_vec, {vertcat(obj.prev_fe.cross_comp_cont_1{end,:})}];
            else
                cross_comp_cont_1_vec = obj.cross_comp_cont_1(:,p.Results.sys).';
                cross_comp_cont_1_vec = [cross_comp_cont_1_vec, obj.prev_fe.cross_comp_cont_1(end,p.Results.sys)];
            end
            sum_cross_comp_cont_1 = sum([cross_comp_cont_1_vec{:}], 2);
        end


        function sum_cross_comp_cont_2 = sumCrossCompCont2(obj, varargin)
            % TODO!!! : DO appropiate cross comp functions
            import casadi.*
            p = inputParser();
            p.FunctionName = 'sumCrossCompCont2';

            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});

            if ismember('sys', p.UsingDefaults)
                cross_comp_cont_2 = obj.cross_comp_cont_2;
                cross_comp_cont_2_vec = arrayfun(@(row) vertcat(cross_comp_cont_2{row, :}), 1:size(cross_comp_cont_2,1), 'UniformOutput', false);
                cross_comp_cont_2_vec = [cross_comp_cont_2_vec, {vertcat(obj.prev_fe.cross_comp_cont_2{end,:})}];
            else
                cross_comp_cont_2_vec = obj.cross_comp_cont_2(:,p.Results.sys).';
                cross_comp_cont_2_vec = [cross_comp_cont_2_vec, obj.prev_fe.cross_comp_cont_2(end,p.Results.sys)];
            end
            sum_cross_comp_cont_2 = sum([cross_comp_cont_2_vec{:}], 2);
        end


        function sum_elastic = sumElastic(obj)
            elastic = obj.elastic;
            sum_elastic = sum(elastic);
        end

        function z = rkStageZ(obj, stage)
            import casadi.*

            % TODO: theta_step/beta
            idx = [[obj.ind_theta{stage, :}],...
                [obj.ind_lam{stage, :}],...
                [obj.ind_mu{stage, :}],...
                [obj.ind_alpha{stage, :}],...
                [obj.ind_lambda_n{stage, :}],...
                [obj.ind_lambda_p{stage, :}],...
                [obj.ind_beta{stage}],...
                [obj.ind_theta_step{stage}],...
                [obj.ind_z{stage}],...
                [obj.ind_lambda_normal{stage}],...
                [obj.ind_lambda_tangent{stage}],...
                [obj.ind_y_gap{stage}],...
                [obj.ind_gamma{stage}],...
                [obj.ind_beta_conic{stage}],...
                [obj.ind_gamma_d{stage}],...
                [obj.ind_beta_d{stage}],...
                [obj.ind_delta_d{stage}],...
                [obj.ind_p_vt{stage}],...
                [obj.ind_n_vt{stage}],...
                [obj.ind_alpha_vt{stage}],...
                ];

            z = obj.w(idx);
        end

        function forwardSimulation(obj, ocp, Uk, s_sot, p_stage)
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            obj.u = Uk;
            % left bondary point
            if settings.dcs_mode == "CLS"
                % do continuity on x and impulse equtions
                X_k0 = obj.w(obj.ind_x_left_bp{1}); % corresponds to post impact t^+
                X_k = obj.prev_fe.x{end}; % corresponds to pre impact t^-
                Q_k0  = X_k0(1:dims.n_q);
                V_k0  = X_k0(dims.n_q+1:end);
                Q_k  = X_k(1:dims.n_q);
                V_k  = X_k(dims.n_q+1:end);
                % junction equations
                obj.addConstraint(Q_k0-Q_k);
                
                Z_impulse_k = [obj.w(obj.ind_Lambda_normal{1}), obj.w(obj.ind_Y_gap{1}),obj.w(obj.ind_P_vn{1}),obj.w(obj.ind_N_vn{1}),...
                               obj.w(obj.ind_Lambda_tangent{1}), obj.w(obj.ind_Gamma_d{1}), obj.w(obj.ind_Beta_d{1}), obj.w(obj.ind_Delta_d{1}), ...
                               obj.w(obj.ind_Gamma{1}), obj.w(obj.ind_Beta_conic{1}), obj.w(obj.ind_P_vt{1}), obj.w(obj.ind_N_vt{1}), obj.w(obj.ind_Alpha_vt{1})];
                obj.addConstraint(model.g_impulse_fun(V_k0,V_k,Z_impulse_k));
            else
                X_k0 = obj.prev_fe.x{end};
            end

            if settings.irk_representation == IrkRepresentation.integral
                X_ki = obj.x;
                Xk_end = settings.D_irk(1) * X_k0;
            elseif settings.irk_representation == IrkRepresentation.differential
                X_ki = {};
                for j = 1:dims.n_s
                    x_temp = X_k0;
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*settings.A_irk(j,r)*obj.v{r};
                    end
                    X_ki = [X_ki {x_temp}];
                end
                X_ki = [X_ki, {obj.x{end}}];
                Xk_end = X_k0;
            elseif settings.irk_representation == IrkRepresentation.differential_lift_x
                X_ki = obj.x;
                Xk_end = X_k0;
                for j = 1:dims.n_s
                    x_temp = X_k0;
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*settings.A_irk(j,r)*obj.v{r};
                    end
                    obj.addConstraint(obj.x{j}-x_temp);
                end
            end

            for j = 1:dims.n_s
                % Multiply by s_sot_k which is 1 if not using speed of time variable
                [fj, qj] = model.f_x_fun(X_ki{j}, obj.rkStageZ(j), Uk, p_stage, model.v_global);
                fj = s_sot*fj;
                qj = s_sot*qj;
                gj = model.g_z_all_fun(X_ki{j}, obj.rkStageZ(j), Uk, p_stage, model.v_global);

                obj.addConstraint(gj);
                if settings.irk_representation == IrkRepresentation.integral
                    xj = settings.C_irk(1, j+1) * X_k0;
                    for r = 1:dims.n_s
                        xj = xj + settings.C_irk(r+1, j+1) * X_ki{r};
                    end
                    Xk_end = Xk_end + settings.D_irk(j+1) * X_ki{j};
                    obj.addConstraint(obj.h * fj - xj);
                    obj.cost = obj.cost + settings.B_irk(j+1) * obj.h * qj;
                    obj.objective = obj.objective + settings.B_irk(j+1) * obj.h * qj;
                else
                    Xk_end = Xk_end + obj.h * settings.b_irk(j) * obj.v{j};
                    obj.addConstraint(fj - obj.v{j});
                    obj.cost = obj.cost + settings.b_irk(j) * obj.h * qj;
                    obj.objective = obj.objective + settings.b_irk(j) * obj.h * qj;
                end
            end

            % nonlinear inequality.
            % TODO: do this cleaner
            if (model.g_path_constraint &&...
                    (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_path_at_fe))
                obj.addConstraint(model.g_path_fun(X_k0,Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
            end
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if model.g_path_constraint && settings.g_path_at_stg
                    obj.addConstraint(model.g_path_fun(X_ki{j},Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
                end
            end

            % end constraints
            if (~settings.right_boundary_point_explicit ||...
                    settings.irk_representation == IrkRepresentation.differential)
                obj.addConstraint(Xk_end - obj.x{end});
            end
            if (~settings.right_boundary_point_explicit &&...
                    settings.use_fesd &&...
                    obj.fe_idx < dims.N_finite_elements(obj.ctrl_idx))

                % TODO verify this.
                obj.addConstraint(model.g_switching_fun(obj.x{end}, obj.rkStageZ(dims.n_s+1), Uk, p_stage));
            end
        end

        function createComplementarityConstraints(obj, sigma_p, s_elastic, p_stage)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            psi_fun = settings.psi_fun;

            % TODO This needs n_comp. that can be calculated apriori so lets do that
            if settings.elasticity_mode == ElasticityMode.ELL_1
                s_elastic = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['s_elastic_' num2str(obj.ctrl_idx) '_' num2str(obj.fe_idx)],...
                    n_comp);
                obj.addVariable(s_elastic,...
                    'elastic',...
                    settings.s_elastic_min*ones(n_comp,1),...
                    settings.s_elastic_max*ones(n_comp,1),...
                    settings.s_elastic_0*ones(n_comp,1));
            end

            if settings.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            else
                sigma = s_elastic;
            end
            
            
            g_path_comp = [];
            % path complementarities
            if (model.g_comp_path_constraint &&...
                (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_path_at_fe))
                pairs = model.g_comp_path_fun(obj.prev_fe.x{end}, obj.u, p_stage, model.v_global);
                expr = apply_psi(pairs, psi_fun, sigma_p);
                g_path_comp = vertcat(g_path_comp, expr);
            end
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if model.g_comp_path_constraint && settings.g_path_at_stg
                    pairs = model.g_comp_path_fun(obj.x{j}, obj.u, p_stage, model.v_global);
                    expr = apply_psi(pairs, psi_fun, sigma_p);
                    g_path_comp = vertcat(g_path_comp, expr);
                end
            end

            g_impulse_comp = [];
            impulse_pairs = [];
            if settings.dcs_mode == DcsMode.CLS
                 % comp condts 
                % Y_gap comp. to Lambda_Normal+P_vn+N_vn;..
                % P_vn comp.to N_vn;
                %  if conic:
                    % Gamma comp. to Beta
                    % if abs: P_vt comp.to N_vt;
                    % if lp:  P_vt comp.to e - Alpha_vt, N_vt comp.to Alpha_vt
                %  if polyhedral
                   % Delta comp. to Lambda_tangent
                   % Gamma comp. to Beta;
                
                Y_gap = obj.w(obj.ind_Y_gap{1});
                Lambda_normal = obj.w(obj.ind_Lambda_normal{1});
                P_vn = obj.w(obj.ind_P_vn{1});
                N_vn = obj.w(obj.ind_N_vn{1});
                P_vt = obj.w(obj.ind_P_vt{1});
                N_vt = obj.w(obj.ind_N_vt{1});
                Gamma = obj.w(obj.ind_Gamma{1});
                Beta_conic = obj.w(obj.ind_Beta_conic{1});
                Lambda_tangent = obj.w(obj.ind_Lambda_tangent{1});
                Beta_d = obj.w(obj.ind_Beta_d{1});
                Alpha_vt = obj.w(obj.ind_Alpha_vt{1});
                Delta_d = obj.w(obj.ind_Delta_d{1});
                Gamma_d = obj.w(obj.ind_Gamma_d{1});
                
                impulse_pairs = vertcat(impulse_pairs, [Y_gap, (Lambda_normal+P_vn+N_vn)]);
                impulse_pairs = vertcat(impulse_pairs, [P_vn, N_vn]);
                if model.friction_exists
                    if settings.friction_model == FrictionModel.Conic
                        g_impulse_comp = vertcat(g_impulse_comp, Gamma*Beta_conic);
                        switch settings.conic_model_switch_handling
                          case ConicModelSwitchHandling.Plain
                            % no extra comps
                          case ConicModelSwitchHandling.Abs
                            impulse_pairs = vertcat(impulse_pairs, [P_vt,N_vt]);
                          case ConicModelSwitchHandling.Lp
                            impulse_pairs = vertcat(impulse_pairs, [P_vt,1-Alpha_vt]);
                            impulse_pairs = vertcat(impulse_pairs, [Alpha_vt,N_vt]);
                        end
                    elseif settings.friction_model == FrictionModel.Polyhedral
                        impulse_pairs = vertcat(impulse_pairs, [Delta_d,Lambda_tangent]);
                        impulse_pairs = vertcat(impulse_pairs, [Gamma_d,Beta_d]);
                    end
                end
                expr = apply_psi(impulse_pairs, psi_fun, sigma_p);
                g_impulse_comp = expr
            end


            % Generate all complementarity pairs
            cross_comp_pairs = cell(dims.n_s, dims.n_s+1, dims.n_sys);
            cross_comp_discont_0 = obj.cross_comp_discont_0;
            cross_comp_cont_0 = obj.cross_comp_cont_0;
            % Complement within FE
            for j=1:dims.n_s
                for jj = 1:dims.n_s
                    for r=1:dims.n_sys
                        cross_comp_pairs{j, jj+1, r} = horzcat(cross_comp_discont_0{j,r},cross_comp_cont_0{jj,r});
                    end
                end
            end
            for j=1:dims.n_s
                for r=1:dims.n_sys
                    cross_comp_pairs{j,1,r} = horzcat(cross_comp_discont_0{j,r}, obj.prev_fe.cross_comp_cont_0{end,r});
                end
            end
            obj.cross_comp_pairs = cross_comp_pairs;

            % apply psi
            g_cross_comp = [];
            if settings.cross_comp_mode == 1
                for j=1:dims.n_s
                    for jj = 1:dims.n_s+1
                        for r=1:dims.n_sys
                            pairs = cross_comp_pairs{j, jj, r};
                            g_cross_comp = vertcat(g_cross_comp, apply_psi(pairs, psi_fun, sigma));
                        end
                    end
                end
            elseif settings.cross_comp_mode == 2
                for j=1:dims.n_s
                    for jj = 1:dims.n_s+1
                        for r=1:dims.n_sys
                            pairs = cross_comp_pairs{j, jj, r};
                            g_cross_comp = vertcat(g_cross_comp, sum(apply_psi(pairs, psi_fun, sigma/size(pairs, 1))));
                        end
                    end
                end
            elseif settings.cross_comp_mode == 3
                for j=1:dims.n_s
                    for r=1:dims.n_sys
                        pairs = cross_comp_pairs(j, :, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/(dims.n_s+1)), pairs, 'uni', false);
                        exprs = sum2([expr_cell{:}]);
                        g_cross_comp = vertcat(g_cross_comp, exprs);
                    end
                end
            elseif settings.cross_comp_mode == 4
                for jj=1:dims.n_s+1
                    for r=1:dims.n_sys
                        pairs = cross_comp_pairs(:, jj, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/(dims.n_s)), pairs, 'uni', false);
                        exprs = sum2([expr_cell{:}]);
                        g_cross_comp = vertcat(g_cross_comp, exprs);
                    end
                end
            elseif settings.cross_comp_mode == 5
                for j=1:dims.n_s
                    for r=1:dims.n_sys
                        pairs = cross_comp_pairs(j, :, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/((dims.n_s+1)*dims.n_theta)), pairs, 'uni', false);
                        exprs = sum1(sum2([expr_cell{:}]));
                        g_cross_comp = vertcat(g_cross_comp, exprs);
                    end
                end
            elseif settings.cross_comp_mode == 6
                for jj=1:dims.n_s+1
                    for r=1:dims.n_sys
                        pairs = cross_comp_pairs(:, jj, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/(dims.n_s*dims.n_theta)), pairs, 'uni', false);
                        exprs = sum1(sum2([expr_cell{:}]));
                        g_cross_comp = vertcat(g_cross_comp, exprs);
                    end
                end
            elseif settings.cross_comp_mode == 7
                for r=1:dims.n_sys
                    pairs = cross_comp_pairs(:, :, r);
                    expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/((dims.n_s+1)*dims.n_s)), pairs, 'uni', false);
                    exprs = sum2([expr_cell{:}]);
                    g_cross_comp = vertcat(g_cross_comp, exprs);
                end
            elseif settings.cross_comp_mode == 8
                for r=1:dims.n_sys
                    pairs = cross_comp_pairs(:, :, r);
                    expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma/((dims.n_s+1)*dims.n_s*dims.n_theta)), pairs, 'uni', false);
                    exprs = sum1(sum2([expr_cell{:}]));
                    g_cross_comp = vertcat(g_cross_comp, exprs);
                end
            elseif settings.cross_comp_mode > 8
                return
            end

            g_comp = vertcat(g_cross_comp, g_path_comp, g_impulse_comp);

            [g_comp_lb, g_comp_ub, g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings);
            
            n_cross_comp = length(g_cross_comp);
            n_path_comp = length(g_path_comp);
            n_comp = n_cross_comp + n_path_comp;

            [g_comp, g_comp_lb, g_comp_ub, cost] = reformulate_complementarities(g_comp, settings.mpcc_mode, sigma_p, s_elastic);

            % Add reformulated constraints
            obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);

            % If We need to add a cost from the reformulation do that as needed;
            if settings.mpcc_mode == MpccMode.ell_1_penalty
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*cost;
                else
                    obj.cost = sigma_p*obj.cost + cost;
                end
            end
        end

        function stepEquilibration(obj, sigma_p, rho_h_p)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            % TODO implement other modes!
            if settings.use_fesd && obj.fe_idx > 1
                nu = obj.nu_vector;
                delta_h_ki = obj.h - obj.prev_fe.h;
                if settings.step_equilibration == StepEquilibrationMode.heuristic_mean
                    h_fe = model.T / (sum(dims.N_finite_elements)); % TODO this may be a bad idea if using different N_fe. may want to issue warning in that case
                    obj.cost = obj.cost + rho_h_p * (obj.h - h_fe).^2;
                elseif settings.step_equilibration ==  StepEquilibrationMode.heuristic_diff
                    obj.cost = obj.cost + rho_h_p * delta_h_ki.^2;
                elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed_scaled
                    obj.cost = obj.cost + rho_h_p * tanh(nu/settings.step_equilibration_sigma) * delta_h_ki.^2;
                elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed
                    obj.cost = obj.cost + rho_h_p * nu * delta_h_ki.^2
                elseif settings.step_equilibration == StepEquilibrationMode.direct
                    obj.addConstraint(nu*delta_h_ki, 0, 0);
                elseif settings.step_equilibration == StepEquilibrationMode.direct_homotopy
                    obj.addConstraint([nu*delta_h_ki-sigma_p;-nu*delta_h_ki-sigma_p],...
                        [-inf;-inf],...
                        [0;0]);
                elseif settings.step_equilibration == StepEquilibrationMode.direct_homotopy_lift
                    obj.addConstraint([obj.nu_lift-nu;obj.nu_lift*delta_h_ki-sigma_p;-obj.nu_lift*delta_h_ki-sigma_p],...
                        [0;-inf;-inf],...
                        [0;0;0]);
                else
                    error("Step equilibration mode not implemented");
                end
            end
        end
    end
end

