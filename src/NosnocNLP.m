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


% TODO cleanup steps:
%      - Create primal variables all at once.
%      - Separate sections into separate functions operating on the `problem` struct/class
%      - time variables should probably not just be lumped into the state, for readability.
%      - remove index in symbolic variable defintions and add instructive
%        names, e.g., Uk -> U,  h_ki -> h_fe, X_ki_stages ->  X_rk_stages
%      - provide instructive names for terminal constraint relaxations
%      - provide more instructive names for cross_comp (match python)

classdef NosnocNLP < NosnocFormulationObject
    properties
        ind_elastic
        ind_c_lift

        % Parameter index variables
        ind_p_x0
        ind_p_global
        ind_p_time_var

        % Problem data
        mpcc
        solver_options

        % original initialization
        w0_original

        % Algorithmic parameters
        sigma_p

        % Algorithmic global variables (time independent)
        s_elastic_inf
        s_elastic_1

        % Parameters
        p
        p0

        % Problem cost function
        augmented_objective_fun

        % Problem objective function
        objective_fun

        % Problem constraint function
        g_fun

        % switch indicator function
        nu_fun

        nabla_J
        nabla_J_fun

        % mapping mpcc indices to nlp indices
        ind_map
    end
    % remaining list of TODOs
    % TODO: cleanup/add properties (in all components)
    % TODO: Create solver object, which will interact with setting parameters.

    properties(Dependent, SetAccess=private, Hidden)
        % Properties generated on the fly.

        % casadi symbolics/expresions for u, sot, and nu
        u
        sot
        nu_vector
        cc_vector

        % Indices for all algebraic vars in the problem
        ind_z_all
        ind_x_all
    end

    methods
        function obj = NosnocNLP(solver_options, mpcc)
            import casadi.*
            obj@NosnocFormulationObject();

            obj.mpcc = mpcc;
            obj.solver_options = solver_options;
            
            sigma_p = define_casadi_symbolic(mpcc.problem_options.casadi_symbolic_mode, 'sigma_p');
            obj.sigma_p = sigma_p;
            obj.p = [sigma_p;mpcc.p];

            obj.augmented_objective = mpcc.augmented_objective;
            obj.objective = mpcc.objective;

            obj.mpcc_to_nlp(sigma_p);
            
            % Process elastic costs
            if solver_options.elasticity_mode == ElasticityMode.ELL_INF
                if solver_options.objective_scaling_direct
                    obj.augmented_objective = obj.augmented_objective + (1/sigma_p)*obj.s_elastic_inf;
                else
                    obj.augmented_objective = sigma_p*obj.augmented_objective + obj.s_elastic_inf;
                end
            end
            if solver_options.elasticity_mode == ElasticityMode.ELL_1
                sum_s_elastic = sum1(obj.w(obj.ind_elastic));
                if solver_options.objective_scaling_direct
                    obj.augmented_objective = obj.augmented_objective + (1/sigma_p)*sum_s_elastic;
                else
                    obj.augmented_objective = sigma_p*obj.augmented_objective + sum_s_elastic;
                end
            end

            % process lifting cost if necessary, TODO: should be parametric
            if solver_options.mpcc_mode == MpccMode.lifting
                c_lift = obj.w(obj.ind_c_lift);

                obj.augmented_objective = obj.augmented_objective + solver_options.lifting_phase1_tau*sum1(c_lift);
            end
            obj.augmented_objective_fun = Function('augmented_objective_fun', {obj.w, obj.p}, {obj.augmented_objective});
            obj.objective_fun = Function('objective_fun', {obj.w, obj.p}, {obj.objective});
            obj.g_fun = Function('g_fun', {obj.w, obj.p}, {obj.g});

            obj.p0 = [solver_options.sigma_0; mpcc.p0];

            obj.w0_original = obj.w0;
           
            % create CasADi function for cost gradient.
            nabla_J = obj.augmented_objective.jacobian(obj.w);
            nabla_J_fun = Function('nabla_J_fun', {obj.w,obj.p},{nabla_J});
            obj.nabla_J = nabla_J;
            obj.nabla_J_fun = nabla_J_fun;
        end

        % TODO deduplicate 
        function mpcc_to_nlp(obj, sigma_p)
            import casadi.*
            mpcc = obj.mpcc;

            % add zeroth finite element variables and constraints (this contains on complementarity pair, for now)
            fe0 = mpcc.fe0;
            X0 = fe0.x{1};
            obj.addPrimalVector(X0,...
                fe0.lbw(fe0.ind_x{1}),...
                fe0.ubw(fe0.ind_x{1}),...
                fe0.w0(fe0.ind_x{1}));
            obj.addConstraint(fe0.g, fe0.lbg, fe0.ubg);

            % Add global vars
            obj.addPrimalVector(mpcc.w(mpcc.ind_v_global),...
                mpcc.lbw(mpcc.ind_v_global),...
                mpcc.ubw(mpcc.ind_v_global),...
                mpcc.w0(mpcc.ind_v_global));

            % Add elastic variable if ell_inf mode
            if obj.solver_options.elasticity_mode == ElasticityMode.ELL_INF
                s_elastic = define_casadi_symbolic(mpcc.problem_options.casadi_symbolic_mode, 's_elastic',1);
                obj.s_elastic_inf = s_elastic;
                if obj.solver_options.elastic_scholtes
                    obj.solver_options.s_elastic_max = inf;
                    obj.addConstraint(s_elastic-obj.sigma_p,-inf,0);
                end
                obj.addVariable(s_elastic,...
                    'elastic',...
                    obj.solver_options.s_elastic_min,...
                    obj.solver_options.s_elastic_max,...
                    obj.solver_options.s_elastic_0);
            else
                s_elastic = [];
            end

            
            if ~mpcc.problem_options.local_speed_of_time_variable && mpcc.problem_options.use_speed_of_time_variables
                sot = mpcc.w(mpcc.ind_sot{1});
                lbsot = mpcc.lbw(mpcc.ind_sot{1});
                ubsot = mpcc.ubw(mpcc.ind_sot{1});
                sot0 = mpcc.w0(mpcc.ind_sot{1});
                obj.addPrimalVector(sot, lbsot, ubsot, sot0);
            end
            for stage=obj.mpcc.stages
                [u, lbu, ubu, u0] = stage.u;
                obj.addPrimalVector(u, lbu, ubu, u0);
                [sot, lbsot, ubsot, sot0] = stage.sot;
                obj.addPrimalVector(sot, lbsot, ubsot, sot0);

                [g_stage, lbg_stage, ubg_stage] = stage.g_stage;
                obj.addConstraint(g_stage, lbg_stage, ubg_stage);
                for fe=stage.stage
                    obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);
                    obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                    if obj.solver_options.mpcc_mode == MpccMode.lifting
                        obj.lift_complementarity_constraints(fe);
                    else
                        obj.relax_complementarity_constraints(fe);
                    end
                end
            end
            % Add s_terminal
            % Add global vars
            obj.addPrimalVector(mpcc.w(mpcc.ind_s_terminal),...
                mpcc.lbw(mpcc.ind_s_terminal),...
                mpcc.ubw(mpcc.ind_s_terminal),...
                mpcc.w0(mpcc.ind_s_terminal));

            % Add t final
            obj.addPrimalVector(mpcc.w(mpcc.ind_t_final),...
                mpcc.lbw(mpcc.ind_t_final),...
                mpcc.ubw(mpcc.ind_t_final),...
                mpcc.w0(mpcc.ind_t_final));

            obj.addConstraint(mpcc.g(mpcc.ind_g_mpcc),...
                mpcc.lbg(mpcc.ind_g_mpcc),...
                mpcc.ubg(mpcc.ind_g_mpcc));

            % Generate index map
            % TODO solver should handle this for clean interface.
            ind_map = 1:length(obj.w);
            ind_map([obj.ind_elastic, obj.ind_c_lift]) = [];
            obj.ind_map = ind_map;
        end

        function relax_complementarity_constraints(obj, component)
            sigma_p = obj.sigma_p;

            psi_fun = obj.solver_options.psi_fun;

            if obj.solver_options.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            elseif obj.solver_options.elasticity_mode == ElasticityMode.ELL_INF
                sigma = obj.s_elastic_inf;
            else
            end
            
            comp_pairs = component.all_comp_pairs;
            n_comp_pairs = size(comp_pairs, 1);

            if obj.solver_options.mpcc_mode == MpccMode.ell_1_penalty
                for ii=1:n_comp_pairs
                    obj.augmented_objective = obj.augmented_objective + comp_pairs(ii,1).*comp_pairs(ii,2);
                end
            else                
                for ii=1:n_comp_pairs
                    if obj.solver_options.elasticity_mode == ElasticityMode.ELL_1
                        n_pairs = length(comp_pairs(ii,1));
                        s_elastic = define_casadi_symbolic(obj.mpcc.problem_options.casadi_symbolic_mode, 's_elastic',n_pairs);
                        sigma = s_elastic;
                        
                        if obj.solver_options.elastic_scholtes
                            obj.solver_options.s_elastic_max = inf;
                            obj.addConstraint(s_elastic-obj.sigma_p,-inf,0);
                        end
                        obj.addVariable(s_elastic,...
                            'elastic',...
                            obj.solver_options.s_elastic_min,...
                            obj.solver_options.s_elastic_max,...
                            obj.solver_options.s_elastic_0);
                    end
                    expr = psi_fun(comp_pairs(ii,1), comp_pairs(ii,2), sigma);
                    [lb, ub, expr] = generate_mpcc_relaxation_bounds(expr, obj.solver_options);
                    obj.addConstraint(expr, lb, ub);
                end
            end
        end

        function lift_complementarity_constraints(obj, component)
            import casadi.*

            a = SX.sym('a');
            b = SX.sym('b');
            c = SX.sym('c');
            s_fun = obj.solver_options.s_fun;

            psi_fun = Function('psi', {a,b,c}, {vertcat(a-s_fun(c), b-s_fun(-c))});

            comp_pairs = component.all_comp_pairs;
            n_comp_pairs = size(comp_pairs, 1);

            for ii=1:n_comp_pairs
                n_pairs = length(comp_pairs(ii,1));
                c_lift = define_casadi_symbolic(obj.mpcc.problem_options.casadi_symbolic_mode, 'c_lift',n_pairs);
              
                obj.addVariable(c_lift,...
                    'c_lift',...
                    -inf,...
                    inf,...
                    0); % TODO a smarter initialization as in Stein2010
                expr = psi_fun(comp_pairs(ii,1), comp_pairs(ii,2), c_lift);
                obj.addConstraint(expr);
            end
        end
        
        function addPrimalVector(obj, symbolic, lb, ub, initial)
            lens = [size(symbolic,1), size(lb,1), size(ub,1), size(initial,1)];
            if ~all(lens == lens(1))
                symbolic
                lb
                ub
                initial
                error("mismatched dims")
            end
            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = vertcat(obj.lbw, lb);
            obj.ubw = vertcat(obj.ubw, ub);
            obj.w0 = vertcat(obj.w0, initial);
        end

        function print(obj,filename)
            if exist('filename')
                delete(filename);
                fileID = fopen(filename, 'w');
            else
                fileID = 1;
            end
            fprintf(fileID, "i\tlbg\t\t ubg\t\t g_expr\n");
            for i = 1:length(obj.lbg)
                expr_str = formattedDisplayText(obj.g(i));
                fprintf(fileID, "%d\t%.2e\t%.2e\t%s\n", i, obj.lbg(i), obj.ubg(i), expr_str);
            end

            fprintf(fileID, "\nw\t\t\tw0\t\tlbw\t\tubw\n");
            for i = 1:length(obj.lbw)
                % keyboard
                expr_str = pad(formattedDisplayText(obj.w(i)), 20);
                lb_str = pad(sprintf('%.2e', obj.lbw(i)), 10);
                fprintf(fileID, "%s\t%.2e\t%s\t%.2e\t\n", expr_str, obj.w0(i), lb_str, obj.ubw(i));
            end

            fprintf(fileID, '\naugmented objective\n');
            fprintf(fileID, strcat(formattedDisplayText(obj.augmented_objective), '\n'));
        end
    end
end
