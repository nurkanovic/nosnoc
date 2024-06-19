function [pss_model] = cls_inelastic(cls_model, opts)
    dims = cls_model.dims;
    pss_model = nosnoc.model.Pss();
    pss_model.dims = dims;

    % transfer all common properties to pss_model
    common_props = properties('nosnoc.model.Base')
    for ii=1:length(common_props)
        prop = common_props{ii};
        pss_model.(prop) = cls_model.(prop);
    end
    
    dims.n_contacts = length(cls_model.f_c);
    % check dimensions of contacts
    if isempty(dims.n_dim_contact)
        warning('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
        dims.n_dim_contact = 2;
    end
    
    % quadrature state
    dims.n_quad = 0;
    f_quad = [];
    if opts.time_freezing_quadrature_state
        % define quadrature state
        L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
        pss_model.lbx = [pss_model.lbx;-inf];
        pss_model.ubx = [pss_model.ubx;inf];
        pss_model.x = [pss_model.x;L];
        pss_model.x0 = [pss_model.x0;0];
        f_quad = pss_model.f_q;
        pss_model.f_q = 0;
        if ~isempty(pss_model.f_q_T)
            pss_model.f_q_T  = pss_model.f_q_T + L;
        else
            pss_model.f_q_T = L;
        end
        dims.n_quad = 1;
    end
    
    % uneven number of states = it is assumed that the clock state is defined.
    t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
    % update lower and upper bounds of lbx and ubx
    pss_model.lbx = [pss_model.lbx;-inf];
    pss_model.ubx = [pss_model.ubx;inf];
    pss_model.x = [pss_model.x;t];
    pss_model.x0 = [pss_model.x0;0];
    
    % normal and tangential velocities
    eps_t = 1e-7;
    v_normal = model.J_normal'*model.v;
    if obj.friction_exists
        if dims.n_dim_contact == 2
            v_tangent = (model.J_tangent'*model.v)';
        else
            v_tangent = model.J_tangent'*model.v;
            v_tangent = reshape(v_tangent,2,dims.n_contacts); % 2 x n_c , the columns are the tangential velocities of the contact points

        end
        v_tangent_norms = [];
        for ii = 1:dims.n_contacts
            v_tangent_norms = [v_tangent_norms;norm(v_tangent(:,ii))];
        end
    else
        v_tangent  = [];
    end
    
    % parameter for auxiliary dynamics
    % TODO add to options
    a_n = 100;
    %% Time-freezing reformulation
    % unconstrained dynamics with clock state
    f_ode = [model.v;...
        model.invM*model.f_v;
        f_quad;
        1];

    % Auxiliary dynamics
    % where to use invM, in every auxiliary dynamics or only at the end
    if inv_M_once
        inv_M_aux = eye(dims.n_q);
        inv_M_ext = blkdiag(zeros(dims.n_q),model.invM,0);
    else
        inv_M_aux = model.invM;
        inv_M_ext = eye(dims.n_x+1);
    end
    f_aux_pos = []; % matrix with all auxiliary tangent dynamics
    f_aux_neg = [];
    % time freezing dynamics
    if opts.stabilizing_q_dynamics
        f_q_dynamics = -opts.kappa_stabilizing_q_dynamics*obj.J_normal*diag(model.f_c);
    else
        f_q_dynamics = zeros(dims.n_q,dims.n_contacts);
    end
    f_aux_normal = [f_q_dynamics;inv_M_aux*obj.J_normal*obj.a_n;zeros(1,dims.n_contacts)];

    if opts.nonsmooth_switching_fun
        pss_model.c = [max_smooth_fun(model.f_c,v_normal,0);v_tangent];    
    else
        if dims.n_dim_contact == 2
            pss_model.c = [model.f_c;v_normal;v_tangent'];
        else
            pss_model.c = [model.f_c;v_normal;v_tangent_norms-eps_t];
        end
    end

    for ii = 1:dims.n_contacts
        if model.friction_exists && model.mu_f(ii)>0
            % auxiliary tangent;
            if dims.n_dim_contact == 2
                v_tangent_ii = model.J_tangent(:,ii)'*model.v;
                f_aux_pos_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(model.J_normal(:,ii)-model.J_tangent(:,ii)*(model.mu_f(ii)))*a_n;0]; % for v>0
                f_aux_neg_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(model.J_normal(:,ii)+model.J_tangent(:,ii)*(model.mu_f(ii)))*a_n;0]; % for v<0
            else
                v_tangent_ii = v_tangent(:,ii);
                f_aux_pos_ii = [f_q_dynamics(:,ii);inv_M_aux*(model.J_normal(:,ii)*a_n-model.J_tangent(:,ii*2-1:ii*2)*model.mu_f(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                f_aux_neg_ii = [f_q_dynamics(:,ii);inv_M_aux*(model.J_normal(:,ii)*a_n+model.J_tangent(:,ii*2-1:ii*2)*model.mu_f(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
            end
            f_aux_pos = [f_aux_pos,f_aux_pos_ii];
            f_aux_neg = [f_aux_neg,f_aux_neg_ii];
        end
    end
    if model.friction_exists
        f_aux = [f_aux_pos,f_aux_neg];
        dims.n_aux = 2*dims.n_contacts;
    else
        f_aux = f_aux_normal;
        dims.n_aux = dims.n_contacts;
    end
    
    % Build logical functions because its easier than hand picking logic.
    alpha = SX.sym('alpha', length(pss_model.c));
    if ~problem_options.nonsmooth_switching_fun
        alpha_q = obj.alpha(1:dims.n_contacts);
        alpha_v_normal = obj.alpha(dims.n_contacts+1:2*dims.n_contacts);
        if obj.friction_exists
            alpha_v_tangent = obj.alpha(2*dims.n_contacts+1:end);
        end
    else
        alpha_qv = obj.alpha(1:dims.n_contacts);
        if obj.friction_exists
            alpha_v_tangent = obj.alpha(dims.n_contacts+1:end);
        end
    end
    alpha_ode = 1;
    alpha_aux = zeros(dims.n_aux)
    for ii = 1:dims.n_contacts
        if problem_options.nonsmooth_switching_fun
            alpha_ode = alpha_ode*alpha_qv(ii);
            if cls_model.friction_exists
                alpha_aux(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                alpha_aux(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii+1)=(1-alpha_qv(ii));
            end
        else
            alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
            if cls_model.friction_exists
                alpha_aux(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(alpha_v_tangent(ii));
                alpha_aux(1+dims.n_contacts+ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
            end
        end
    end
    indicator = Function('alpha_indicator' {alpha}, {alpha_ode, alpha_aux});
    
    S = dec2bin(0:2^((2+cls_model.friction_exists)*dims.n_contacts)-1) - '0';
    

    
    pss_model.F = [f_ode (inv_M_ext*f_aux)];

    
    pss_model.S = ones(size(obj.F,2),length(obj.c)); % dummy value to pass error checks
                                                     % number of auxiliary dynamics modes
    
    pss_model.dims = dims;
end
