classdef Pds < nosnoc.model.Base
    properties
        f_x
        c
        E
    end

    methods
        function obj = Pds()
            obj = obj@nosnoc.model.Base();
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            verify_and_backfill@nosnoc.model.Base(obj,opts);

            dims = obj.dims;

            if size(obj.f_x,1) ~= dims.n_x
                error("nosnoc: f_x has incorrect dimension. It must have the same dimension as x.")
            end

            dims.n_c = size(obj.c,1);

            if ~isempty(obj.E)
                if size(obj.E, 1) ~= dims.n_x
                    error("nosnoc: Projection matrix E must be an n_x by n_x matrix where n_x is the number of functions defining the set C.")
                end
                if size(obj.E, 2) ~= dims.n_x
                    error("nosnoc: Projection matrix E must be an n_x by n_x matrix where n_x is the number of functions defining the set C.")
                end
            else
                obj.E = eye(dims.n_c);
            end

            obj.dims = dims;
        end
    end
end
