function [pss_model] = cls_elastic(cls_model, opts)
    dims = cls_model.dims;
    pss_model = nosnoc.model.Pss();
    pss_model.dims = dims;

    % transfer all common properties to pss_model
    common_props = properties('nosnoc.model.Base');
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
        L = define_casadi_symbolic(opts.casadi_symbolic_mode,'L',1);
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
    t = define_casadi_symbolic(opts.casadi_symbolic_mode,'t',1);
    % update lower and upper bounds of lbx and ubx
    pss_model.lbx = [pss_model.lbx;-inf];
    pss_model.ubx = [pss_model.ubx;inf];
    pss_model.x = [pss_model.x;t];
    pss_model.x0 = [pss_model.x0;0];
    
    % normal and tangential velocities
    v_normal = cls_model.J_normal'*cls_model.v;

    % parameter for auxiliary dynamics
    % TODO add to options
    a_n = opts.a_n;
    %% Time-freezing reformulation
    % elastic
    % TODO: add to options
    k_aux = opts.k_aux;
    temp1 = 2*abs(log(cls_model.e));
    temp2 = k_aux/(pi^2+log(cls_model.e)^2);
    c_aux = temp1/sqrt(temp2);
    K = [0 1;-k_aux -c_aux];
    N  = [cls_model.J_normal zeros(dims.n_q,1);...
        zeros(dims.n_q,1) cls_model.invM*cls_model.J_normal];
    f_aux_n1 = N*K*N'*[cls_model.q;cls_model.v];
    f_aux_n1 = [f_aux_n1;zeros(dims.n_quad+1,1)];
    f_ode = [cls_model.v;cls_model.invM*cls_model.f_v;f_quad;1];
    % updated with clock state
    pss_model.F = [f_ode, f_aux_n1];
    pss_model.S = [1; -1];
    pss_model.c = cls_model.f_c;
    dims.n_aux = 1;

    pss_model.dims = dims;
end
