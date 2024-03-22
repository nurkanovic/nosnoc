classdef Model
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

        % TODO: maybe a separate OCP class common to all model types?
        % Objective
        f_q
        f_q_T

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
            
        end
        
        function verify_and_backfill(obj)
            import casadi.*

            dims = obj.dims

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
                obj.u = define_casadi_symbolic(problem_options.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if problem_options.print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
            end
            
            % if ~isfield(obj.data, 'g_path')
            %     obj.data.g_path = [];
            % end
            % if ~isfield(obj.data, 'lbg_path')
            %     obj.data.lbg_path = zeros(size(obj.data.g_path));
            % end
            % if ~isfield(obj.data, 'ubg_path')
            %     obj.data.ubg_path = zeros(size(obj.data.g_path));
            % end
            % if ~isfield(obj.data, 'partial_proj_matrix')
            %     obj.data.partial_proj_matrix = eye(length(obj.data.x));
            % end
        end
    end
end
