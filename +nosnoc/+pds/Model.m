classdef Model < handle
    properties (Access=public)
        % Differential state
        x
        lbx
        ubx
        x0

        % Controls
        u
        lbu
        ubu
        u0

        % Algebraics
        lambda

        % dynamics
        f_x

        % Feasible set
        c

        % ePDS WARNING EXPERIMENTAL
        E

        % user algebraics
        z
        z0
        lbz
        ubz
        g_z % uzer algebraic constraints

        % global variables
        v_global
        v0_global
        lbv_global
        ubv_global
        
        % global parameters
        p_global
        p_global_val

        % time varying parameters
        p_time_var
        p_time_var_val
        % params
        p_time_var_stages
        % All params
        p

        % TODO: maybe a separate OCP class common to all model types?
        % Objective
        f_q
        f_q_T

        % least squares
        lsq_x
        x_ref
        f_lsq_x
        x_ref_val
        lsq_u
        u_ref
        f_lsq_u
        u_ref_val
        lsq_T
        x_ref_end
        f_lsq_T
        x_ref_end_val

        % Path constraints
        g_path
        lbg_path
        ubg_path

        % Terminal constraints
        g_terminal
        lbg_terminal
        ubg_terminal

        % Dimensions
        dims
    end

    methods (Access=public)
        function obj = Model()
            dims = struct;
        end

        function generate_equations(obj)
            
        end

        function generate_variables(obj)
            import casadi.*
            obj.lambda = SX.sym('lambda', obj.dims.n_c);
        end
        
        function verify_and_backfill(obj, opts)
            import casadi.*

            dims = obj.dims;

            %% Check differential state
            if size(obj.x, 1) ~= 0
                dims.n_x = length(obj.x);
                % check  lbx
                if size(obj.lbx, 1) ~= 0
                    if length(obj.lbx) ~= dims.n_x
                        error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    obj.lbx = -inf*ones(dims.n_x,1);
                end
                % check ubx
                if size(obj.ubx, 1) ~= 0
                    if length(obj.ubx) ~= dims.n_x
                        error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    obj.ubx = inf*ones(dims.n_x,1);
                end

                % check x0
                if size(obj.x0, 1) ~= 0
                    if length(obj.x0) ~= dims.n_x
                        error('nosnoc: The vector x0, for the initial value of x has the wrong size.')
                    end
                else
                    obj.x0 = zeros(dims.n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end

            %% Check c
            if size(obj.c,1) ~= 0
                dims.n_c = size(obj.c,1);
            else
                error('nosnoc: Please provide the definition of the active set C.')
            end

            %% Check is u provided
            if size(obj.u, 1) ~= 0
                dims.n_u = length(obj.u);
                % check  lbu
                if size(obj.lbu, 1) ~= 0
                    if length(obj.lbu) ~= dims.n_u
                        error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    obj.lbu = -inf*ones(dims.n_u,1);
                end
                % check ubu
                if size(obj.ubu, 1) ~= 0
                    if length(obj.ubu) ~= dims.n_u
                        error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    obj.ubu = inf*ones(dims.n_u,1);
                end
                % check u0
                if size(obj.u0, 1) ~= 0
                    if length(obj.u0) ~= dims.n_u
                        error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    obj.u0 = 0*ones(dims.n_u,1);
                end
            else
                obj.u = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if opts.print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
            end
            %% Parameters (time variable and that do not change with time)
            if size(obj.p_global, 1) ~= 0
                dims.n_p_global = size(obj.p_global,1);
                if size(obj.p_global_val, 1) ~= 0
                    if size(obj.p_global_val,1) ~= dims.n_p_global
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(dims.n_p_global,1);
                end
            else
                dims.n_p_global = 0;
                obj.p_global = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, opts.N_stages]
                        error('nosnoc: User provided p_time_var_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, opts.N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:opts.N_stages
                    var_full = define_casadi_symbolic(opts.casadi_symbolic_mode, ['p_time_var_' num2str(ii)], dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(opts.casadi_symbolic_mode,'', [0, opts.N_stages]);
                obj.p_time_var_val = double.empty(0,opts.N_stages);
            end

            obj.p = vertcat(obj.p_global,obj.p_time_var);

            %% Stage and terminal costs check
            if size(obj.f_q, 1) == 0
                obj.f_q = 0;
            end

            if size(obj.f_q_T, 1) == 0
                obj.f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if size(obj.lsq_x, 1) ~= 0
                if length(obj.lsq_x)<3
                    error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_x{2},1)~=size(obj.lsq_x{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(obj.lsq_x{1},1)~=size(obj.lsq_x{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
                end

                n_x_ref_rows = size(obj.lsq_x{2},1);
                n_x_ref_cols = size(obj.lsq_x{2},2);
                if n_x_ref_cols == opts.N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    obj.lsq_x{2} = repmat(obj.lsq_x{2},1,opts.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,opts.N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropriate size.')
                end
                obj.x_ref_val = obj.lsq_x{2};
                obj.x_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                obj.f_lsq_x = (obj.lsq_x{1}-obj.x_ref)'*obj.lsq_x{3}*(obj.lsq_x{1}-obj.x_ref);
            else
                obj.x_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref',1);
                obj.f_lsq_x = 0;
                obj.x_ref_val = zeros(1,opts.N_stages);
            end

            % least square terms for control inputs
            if size(obj.lsq_u, 1) ~= 0
                if length(obj.lsq_u)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_u{2},1)~=size(obj.lsq_u{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
                end
                if size(obj.lsq_u{1},1)~=size(obj.lsq_u{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
                end
                n_u_ref_rows = size(obj.lsq_u{2},1);
                n_u_ref_cols = size(obj.lsq_u{2},2);
                if n_u_ref_cols == opts.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    obj.lsq_u{2} = repmat(obj.lsq_u{2},1,opts.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,opts.N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropriate size.')
                end
                obj.u_ref_val = obj.lsq_u{2};
                obj.u_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                obj.f_lsq_u = (obj.lsq_u{1}-obj.u_ref)'*obj.lsq_u{3}*(obj.lsq_u{1}-obj.u_ref);
            else
                obj.u_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'u_ref',1);
                obj.f_lsq_u = 0;
                obj.u_ref_val = zeros(1,opts.N_stages);
            end


            % least square terms for control inputs
            if size(obj.lsq_T, 1) ~= 0
                % sanity chkecs on the input
                if length(obj.lsq_T)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference or the weight matrix are missing.')
                end
                if size(obj.lsq_T{2},1)~=size(obj.lsq_T{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
                end
                if size(obj.lsq_T{1},1)~=size(obj.lsq_T{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the terminal cost do not match.')
                end

                n_x_T_rows = size(obj.lsq_T{2},1);
                n_x_T_cols = size(obj.lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(obj.lsq_T{1}));
                    error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropriate size.')
                end
                obj.x_ref_end_val = obj.lsq_T{2};
                obj.x_ref_end = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                obj.f_lsq_T = (obj.lsq_T{1}-obj.x_ref_end)'*obj.lsq_T{3}*(obj.lsq_T{1}-obj.x_ref_end);
            else
                obj.x_ref_end  = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref_end',1);
                obj.f_lsq_T = 0;
                obj.x_ref_end_val = 0;
            end

            %% Down projection
            if ~isempty(obj.E)
                if size(obj.E) ~= [dims.n_x,dims.n_x]
                    error('nosnoc: Projection matrix E must be of size n_x by n_x')
                end
            else
                obj.E = eye(dims.n_x);
            end
            obj.dims = dims;
        end
    end
end
