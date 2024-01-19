% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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

% TODO prehaps it is necessary to have a separate NLP relaxation and direct mpcc solver class.
classdef NosnocSolver < handle
    properties
        mpcc     % Nosnoc Mpcc
        solver_options  % Nosnoc Options
        nlp   % Nosnoc reformulated NLP

        solver % CasADi solver function

        plugin % Plugin object handling specific solver types

        p_val % parameter values (homotopy parameters + problem parameters)
    end

    methods
        function obj = NosnocSolver(mpcc, solver_options)
            import casadi.*
            tic

            solver_options.preprocess();
            
            obj.mpcc = mpcc;
            obj.solver_options = solver_options;

            switch solver_options.solver
                case 'ipopt'
                    obj.plugin = NosnocIpopt();
                case 'snopt'
                    obj.plugin = NosnocSNOPT();
                case 'worhp'
                    obj.plugin = NosnocWORHP();
                case 'uno'
                    obj.plugin = NosnocUNO();
            end
            
            if solver_options.solver_type == 'RELAXATION_HOMOTOPY'
                obj = obj.construct_relaxation_homotopy();
            elseif solver_options.solver_type == 'DIRECT'
                obj = obj.construct_direct();
            end            
        end

        function set(obj, type, val)
            nlp = obj.nlp;
            mpcc = obj.mpcc;
            if strcmp(type, 'w0')
                if length(obj.nlp.w0) == length(val)
                    obj.nlp.w0 = val;
                else
                    error("nosnoc: if initializing w0 all at once you need to provide a vector of corresponding size.")
                end
            else
                % This line uses the index sets collected during the creation of the NosnocProblem and automatically gets
                % the one of the form 'ind_<type>'. This makes this set generic for all variable types.
                ind = obj.mpcc.(strcat('ind_', type));
                if iscell(ind)
                    flat_ind = sort([ind{:}]);
                else
                    flat_ind = ind;
                end
                % TODO: Interpolation

                if iscell(val)
                    % If the passed value is an N_stage by 1 cell array we assume this initialization is done stage wise
                    if ismatrix(val) && size(val, 1) == mpcc.problem_options.N_stages && size(val,2) == 1
                        for ii=1:mpcc.problem_options.N_stages
                            % All variables of each stage are set to the same value
                            for v=ind(ii,:,:)
                                % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, etc.)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.nlp.w0(nlp.ind_map(v{1})) = val{ii};
                                end
                            end
                        end
                    % Otherwise if we have an initialization of the form N_stages-by-N_fe we do the same but finite-element-wise
                    elseif ismatrix(val) && size(val, 1) == mpcc.problem_options.N_stages && size(val, 2) == mpcc.problem_options.N_finite_elements
                        for ii=1:mpcc.problem_options.N_stages
                            for jj=1:mpcc.problem_options.N_finite_elements
                                % All variables of each finite element are set to the same value
                                for v=ind(ii,jj,:)
                                    % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, cls, etc.)
                                    if ~isempty(v) && length(v{1}) == length(val{ii,jj})
                                        obj.nlp.w0(nlp.ind_map(v{1})) = val{ii,jj};
                                    end
                                end
                            end
                        end
                    else
                        error("nosnoc: Initialization either needs to be Stagewise or Finite Element wise, i.e. provide either an N_stage-by-1 or N_stage-by-N_fe cell array");
                    end
                    % Otherwise we assume that we are initializing via a flat array and we simply check for the same length
                else
                    if ismatrix(val) && length(val) == length(flat_ind)
                        obj.nlp.w0(nlp.ind_map(flat_ind)) = val;
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end

        function [results,stats] = solve(obj)
            if obj.solver_options.solver_type == 'RELAXATION_HOMOTOPY'
                if isstruct(obj.nlp)
                    [results,stats] = obj.solve_generic_homotopy();
                else
                    [results,stats] = obj.solve_relaxation_homotopy();
                end
            elseif obj.solver_options.solver_type == 'DIRECT'
                [results,stats] = obj.solve_direct();
            end
        end

        function violation = compute_constraint_violation(obj, w, g_fun)
            % TODO also handle direct solver
            nlp = obj.nlp;
            ubw_violation = max(max(w - nlp.ubw), 0);
            lbw_violation = max(max(nlp.lbw - w), 0);
            g_val = full(g_fun(w, obj.p_val));
            ubg_violation = max(max(g_val - nlp.ubg), 0);
            lbg_violation = max(max(nlp.lbg - g_val), 0);
            violation = max([lbg_violation, ubg_violation, lbw_violation, ubw_violation]);
        end

        function met = complementarity_tol_met(obj, stats)
            last_stats = stats.solver_stats(end);
            met = 0;
            if abs(stats.complementarity_stats(end)) < 10 * obj.solver_options.comp_tol
                met = 1;
            end
        end

        function print_infeasibility(obj, results)
            if obj.solver_options.print_details_if_infeasible
                print_problem_details(results,obj.model,obj.problem, []);
            end
            if obj.solver_options.pause_homotopy_solver_if_infeasible
                keyboard
            end
        end

        function compute_initial_parameters(obj)
            model = obj.mpcc.model;
            solver_options = obj.solver_options;

            x0 = obj.nlp.w0(1:model.dims.n_x);
            lambda00 = [];
            gamma_00 = [];
            p_vt_00 = [];
            n_vt_00  = [];
            gamma_d00 = [];
            delta_d00 = [];
            y_gap00 = [];
            switch obj.mpcc.problem_options.dcs_mode
              case 'Stewart'
                lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
              case 'Step'
                lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
              case 'CLS'
                % TODO: reconsider this if 0th element has an impulse
                y_gap00 = max(0, model.f_c_fun(x0));
                if model.friction_exists
                    switch obj.mpcc.problem_options.friction_model
                      case 'Polyhedral'
                        v0 = x0(model.dims.n_q+1:end);
                        D_tangent_0 = model.D_tangent_fun(x0);
                        v_t0 = D_tangent_0'*v0;
                        for ii = 1:model.dims.n_contacts
                            ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                            gamma_d00 = [gamma_d00;norm(v_t0(ind_temp))/model.dims.n_t];
                            delta_d00 = [delta_d00;D_tangent_0(:,ind_temp)'*v0+gamma_d00(ii)];
                        end
                      case 'Conic'
                        v0 = x0(model.dims.n_q+1:end);
                        v_t0 = model.J_tangent_fun(x0)'*v0;
                        for ii = 1:model.dims.n_contacts
                            ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                            v_ti0 = v0(ind_temp);
                            gamma_00 = [gamma_00;norm(v_ti0)];
                            switch obj.mpcc.problem_options.conic_model_switch_handling
                              case 'Plain'
                                % no extra vars
                              case {'Abs','Lp'}
                                p_vt_00 = [p_vt_00; max(v_ti0,0)];
                                n_vt_00 = [n_vt_00; max(-v_ti0,0)];
                            end
                        end
                    end
                end
            end
            obj.p_val = [obj.nlp.p0(:);x0(:);lambda00(:);y_gap00(:);gamma_00(:);gamma_d00(:);delta_d00(:);p_vt_00(:);n_vt_00(:)];
        end

        function print_nlp_iter_info(obj, stats)
            solver_stats = stats.solver_stats(end);
            ii = size(stats.solver_stats, 2);

            if strcmp(obj.solver_options.solver, 'ipopt')
                if isfield(solver_stats, 'iterations') && ~isempty(solver_stats.iterations)
                    inf_pr = solver_stats.iterations.inf_pr(end);
                    inf_du = solver_stats.iterations.inf_du(end);
                    objective = solver_stats.iterations.obj(end);
                else
                    inf_pr = nan;
                    inf_du = nan;
                    objective = nan;
                end
                fprintf('%d\t%6.2e\t %6.2e\t %6.2e\t %6.2e \t %6.2e \t %6.3f \t %d \t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                    objective, stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status);
            elseif strcmp(obj.solver_options.solver, 'snopt')
                % TODO: Findout snopt prim du inf log!
                inf_pr = nan;
                inf_du = nan;
                fprintf('%d\t%6.2e\t %6.2e\t %6.2e \t %6.2e \t %6.2e \t %6.3f \t %s \t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                    stats.objective(end), stats.cpu_time(end), solver_stats.secondary_return_status, solver_stats.return_status);
                % TODO: add missing log information
            end
        end

        function print_iterate(obj, iterate, only_violations, filename)
            if exist('filename')
                delete(filename);
                fileID = fopen(filename, 'w');
            else
                fileID = 1;
            end
            if ~exist('only_violations')
                only_violations = 0;
            end
            fileID = 1;
            fprintf(fileID, "\nw\t\t\tlbw\t\tubw\titerate\n");
            for i = 1:length(obj.nlp.lbw)
                if ~only_violations || (iterate(i) < obj.nlp.lbw(i) || iterate(i) > obj.nlp.ubw(i))
                    expr_str = pad(formattedDisplayText(obj.nlp.w(i)), 20);
                    lb_str = pad(sprintf('%.2e', obj.nlp.lbw(i)), 10);
                    ub_str = pad(sprintf('%.2e', obj.nlp.ubw(i)), 10);
                    iterate_str = pad(sprintf('%.2e', iterate(i)), 10);
                    fprintf(fileID, "%s\t%s\t%s\t%s\n", expr_str, lb_str, ub_str, iterate_str);
                end
            end

            % constraints
            g_val = full(obj.nlp.g_fun(iterate, obj.p_val));
            fprintf(fileID, "\ni\tlbg\t\t ubg\t\t g_val\t\tg_expr\n");
            for i = 1:length(obj.nlp.lbg)
                if ~only_violations || (g_val(i) < obj.nlp.lbg(i) || g_val(i) > obj.nlp.ubg(i))
                    expr_str = formattedDisplayText(obj.nlp.g(i));
                    lb_str = pad(sprintf('%.2e', obj.nlp.lbg(i)), 12);
                    ub_str = pad(sprintf('%.2e', obj.nlp.ubg(i)), 12);
                    fprintf(fileID, "%d\t%s\t%s\t%.2e\t%s\n", i, lb_str, ub_str, g_val(i), expr_str);
                end
            end
        end

        function [polished_w, res_out, stat_type, n_biactive] = calculate_stationarity(obj, results, exitfast, lifted)
            import casadi.*
            stat_type = "?";
            nlp = obj.nlp;
            mpcc = obj.mpcc;

            % This method generates a TNLP (tightened NLP) using the results of a homotopy solver.
            % In order to do so we first identify the active set of the complementarity constraints which involves identifying which
            % branch of the complementarity each element of G(w) and H(w) are, or whether the solution is biactive at that point.
            % For a more in depth explanation see the PhD theses of Armin Nurkanovic or Alexandra Schwartz
            
            w_orig = results.x(nlp.ind_map);
            lam_x = full(results.lam_x(nlp.ind_map));
            lam_g = full(results.lam_g);
            lam_g(nlp.ind_g_comp) = []; % TODO: also elastic scholtes constraint

            G = mpcc.G_fun(mpcc.w, mpcc.p);
            H = mpcc.H_fun(mpcc.w, mpcc.p);
            G_old = full(mpcc.G_fun(w_orig, obj.p_val(2:end)));
            H_old = full(mpcc.H_fun(w_orig, obj.p_val(2:end)));
            if obj.solver_options.elasticity_mode == ElasticityMode.NONE
                a_tol = 1*sqrt(obj.solver_options.comp_tol);
            else
                a_tol = 1*sqrt(full(results.x(obj.nlp.ind_elastic)));
            end
            %a_tol = obj.solver_options.comp_tol^2;
            %a_tol = 1e-5;
            ind_00 = G_old<a_tol & H_old<a_tol;
            n_biactive = sum(ind_00)
            if exitfast && n_biactive == 0
                polished_w = [];
                res_out = []
                stat_type = "S";
                disp("Converged to S-stationary point")
                return
            end
            mindists = max(G_old, H_old);
            [~, min_idx] = sort(mindists);

            if lifted
                w = mpcc.w; lbw_orig = mpcc.lbw; ubw_orig = mpcc.ubw;
                f = mpcc.augmented_objective;

                ind_g = setdiff(1:length(nlp.g),nlp.ind_g_comp);
                g = nlp.g(ind_g); lbg = nlp.lbg(ind_g); ubg = nlp.ubg(ind_g);
                p = mpcc.p;

                % Debug
                jac_f = Function('jac_f', {mpcc.w, mpcc.p}, {mpcc.augmented_objective.jacobian(mpcc.w)});
                jac_g = Function('jac_g', {mpcc.w, mpcc.p}, {g.jacobian(mpcc.w)});
                jac_G = Function('jac_G', {mpcc.w, mpcc.p}, {G.jacobian(mpcc.w)});
                jac_H = Function('jac_H', {mpcc.w, mpcc.p}, {H.jacobian(mpcc.w)});

                nabla_f = jac_f(w_orig, obj.p_val(2:end));
                nabla_g = jac_g(w_orig, obj.p_val(2:end));
                nabla_G = jac_G(w_orig, obj.p_val(2:end));
                nabla_H = jac_H(w_orig, obj.p_val(2:end));
                
                lift_G = SX.sym('lift_G', length(G));
                lift_H = SX.sym('lift_H', length(H));
                g_lift = vertcat(lift_G - G, lift_H - H);
                
                w = vertcat(w, lift_G, lift_H);

                lam_x_aug = vertcat(lam_x, zeros(size(G)), zeros(size(H)));
                
                g = vertcat(g,g_lift);
                lbg = vertcat(lbg,zeros(size(g_lift)));
                ubg = vertcat(ubg,zeros(size(g_lift)));
                lam_g_aug = vertcat(lam_g, zeros(size(g_lift)));

                ind_mpcc = 1:length(lbw_orig);
                ind_G = (1:length(lift_G))+length(lbw_orig);
                ind_H = (1:length(lift_H))+length(lbw_orig)+length(lift_G);

                converged = false;
                n_max_biactive = n_biactive;
                n_biactive = 0;
                while ~converged && n_biactive <=n_max_biactive
                    idx_00 = min_idx(1:n_biactive)
                    ind_00 = false(length(G),1);
                    ind_00(idx_00) = true;
                    
                    ind_0p = G_old<a_tol & ~ind_00;
                    ind_p0 = H_old<a_tol & ~ind_00;

                    lblift_G = -inf*ones(size(lift_G));
                    ublift_G = inf*ones(size(lift_G));
                    lblift_H = -inf*ones(size(lift_G));
                    ublift_H = inf*ones(size(lift_G));
                    %lblift_G = zeros(size(lift_G));
                    %ublift_G = inf*ones(size(lift_G));
                    %lblift_H = zeros(size(lift_G));
                    %ublift_H = inf*ones(size(lift_G));
                    ublift_G(find(ind_00 | ind_0p)) = 0; 
                    ublift_H(find(ind_00 | ind_p0)) = 0;

                    G_init = G_old;
                    G_init(ind_00 | ind_0p) = 0;
                    H_init = H_old;
                    H_init(ind_00 | ind_p0) = 0;

                    lbw = vertcat(lbw_orig, lblift_G, lblift_H);
                    ubw = vertcat(ubw_orig, ublift_G, ublift_H);
                    w_init = vertcat(w_orig, G_init, H_init);

                    
                    casadi_nlp = struct('f', f , 'x', w, 'g', g, 'p', p);
                    opts_casadi_nlp.ipopt.max_iter = 5000;
                    opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
                    opts_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
                    opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_bound_frac = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.tol = obj.solver_options.comp_tol;
                    opts_casadi_nlp.ipopt.acceptable_tol = sqrt(obj.solver_options.comp_tol);
                    opts_casadi_nlp.ipopt.acceptable_dual_inf_tol = sqrt(obj.solver_options.comp_tol);
                    opts_casadi_nlp.ipopt.fixed_variable_treatment = 'relax_bounds';
                    opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
                    opts_casadi_nlp.ipopt.bound_relax_factor = 1e-12;
                    opts_casadi_nlp.ipopt.linear_solver = obj.solver_options.opts_casadi_nlp.ipopt.linear_solver;
                    opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
                    opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
                    opts_casadi_nlp.ipopt.nlp_scaling_method = 'none';
                    default_tol = 1e-7;
                    opts_casadi_nlp.ipopt.tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.resto_failure_feasibility_threshold = 0;
                    opts_casadi_nlp.ipopt.print_level = 5;
                    opts_casadi_nlp.ipopt.sb = 'yes';
                    
                    tnlp_solver = nlpsol('tnlp', 'ipopt', casadi_nlp, opts_casadi_nlp);
                    tnlp_results = tnlp_solver('x0', w_init,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'lam_g0', lam_g_aug,...
                        'lam_x0', lam_x_aug,...
                        'p',obj.p_val(2:end));
                    res_out = tnlp_results;
                    if tnlp_solver.stats.success
                        converged = true;
                    else
                        n_biactive = n_biactive + 1;
                    end 
                    %converged = true;
                end
                w_star = full(tnlp_results.x);
                w_star_orig = w_star(ind_mpcc);
                max(abs(w_init - w_star))
                [~,sort_idx] = sort(full(abs(w_init - w_star))); % sorted by difference from original stationary point

                nu = -full(tnlp_results.lam_x(ind_G));
                xi = -full(tnlp_results.lam_x(ind_H));

                % Debug
                jac_f = Function('jac_f', {w, mpcc.p}, {f.jacobian(w)});
                jac_g = Function('jac_g', {w, mpcc.p}, {g.jacobian(w)});
                jac_G = Function('jac_G', {w}, {lift_G.jacobian(w)});
                jac_H = Function('jac_H', {w}, {lift_H.jacobian(w)});

                nabla_f = full(jac_f(w_star, obj.p_val(2:end)));
                nabla_g = full(jac_g(w_star, obj.p_val(2:end)));
                nabla_G = full(jac_G(w_star));
                nabla_H = full(jac_H(w_star));
                
                % Index sets
                I_G = find(ind_0p+ind_00);
                I_H = find(ind_p0+ind_00);

                G_new = full(mpcc.G_fun(w_star_orig, obj.p_val(2:end)));
                H_new = full(mpcc.H_fun(w_star_orig, obj.p_val(2:end)));

                g_tnlp = full(tnlp_results.g);
                % multipliers (nu, xi already calculated above)
                lam_x_star = full(tnlp_results.lam_x);
                lam_x_star(ind_G) = 0; % Capturing multipliers for non-complementarity 
                lam_x_star(ind_H) = 0; % box constraints
                lam_g_star = tnlp_results.lam_g;
                type_tol = a_tol^2;
            else
                w = mpcc.w; lbw_orig = mpcc.lbw; ubw_orig = mpcc.ubw;
                f = mpcc.augmented_objective;

                ind_g = setdiff(1:length(nlp.g),nlp.ind_g_comp);
                g = nlp.g(ind_g); lbg = nlp.lbg(ind_g); ubg = nlp.ubg(ind_g);
                p = mpcc.p;

                lam_x_aug = lam_x;
                
                g = vertcat(g,10000*G,10000*H);
                lbg = vertcat(lbg,zeros(size(G)),zeros(size(H)));
                ubg = vertcat(ubg,zeros(size(G)),zeros(size(H)));
                lam_g_aug = vertcat(lam_g, zeros(size(G)), zeros(size(H)));

                jac_g = Function('jac_g', {mpcc.w, mpcc.p}, {g.jacobian(mpcc.w)});

                ind_mpcc = 1:length(lam_g);
                ind_G = (1:length(G))+length(lam_g);
                ind_H = (1:length(H))+length(lam_g)+length(G);

                converged = false;
                n_max_biactive = n_biactive;
                n_biactive = 0;
                while ~converged && n_biactive <=n_max_biactive
                    idx_00 = min_idx(1:n_biactive)
                    ind_00 = false(length(G),1);
                    ind_00(idx_00) = true;
                    
                    ind_0p = G_old<H_old & ~ind_00;
                    ind_p0 = H_old<G_old & ~ind_00;

                    lbG = 0*ones(size(G));
                    ubG = inf*ones(size(G));
                    lbH = 0*ones(size(G));
                    ubH = inf*ones(size(G));
                    %lbG = zeros(size(G));
                    %ubG = inf*ones(size(G));
                    %lbH = zeros(size(G));
                    %ubH = inf*ones(size(G));
                    ubG(find(ind_00 | ind_0p)) = 1e-12; 
                    ubH(find(ind_00 | ind_p0)) = 1e-12;
                    %g(find(ind_00)) = 10000*g(find(ind_00));
                    lbg(ind_G) = lbG;
                    lbg(ind_H) = lbH;
                    ubg(ind_G) = ubG;
                    ubg(ind_H) = ubH;

                    nabla_g = jac_g(w_orig, obj.p_val(2:end));

                    lbw = lbw_orig;
                    ubw = ubw_orig;
                    w_init = w_orig;

                    
                    casadi_nlp = struct('f', f , 'x', w, 'g', g, 'p', p);
                    opts_casadi_nlp.ipopt.max_iter = 5000;
                    opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
                    %opts_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
                    opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_bound_frac = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.tol = obj.solver_options.comp_tol;
                    opts_casadi_nlp.ipopt.acceptable_tol = sqrt(obj.solver_options.comp_tol);
                    opts_casadi_nlp.ipopt.acceptable_dual_inf_tol = sqrt(obj.solver_options.comp_tol);
                    opts_casadi_nlp.ipopt.fixed_variable_treatment = 'relax_bounds';
                    opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
                    opts_casadi_nlp.ipopt.bound_relax_factor = 1e-12;
                    opts_casadi_nlp.ipopt.linear_solver = obj.solver_options.opts_casadi_nlp.ipopt.linear_solver;
                    opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
                    opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
                    %opts_casadi_nlp.ipopt.nlp_scaling_method = 'equilibration-based';
                    default_tol = 1e-4;
                    opts_casadi_nlp.ipopt.tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
                    %opts_casadi_nlp.ipopt.resto_failure_feasibility_threshold = 0;
                    opts_casadi_nlp.ipopt.print_level = 5;
                    opts_casadi_nlp.ipopt.sb = 'yes';
                    
                    tnlp_solver = nlpsol('tnlp', 'ipopt', casadi_nlp, opts_casadi_nlp);
                    tnlp_results = tnlp_solver('x0', w_init,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'lam_g0', lam_g_aug,...
                        'lam_x0', lam_x_aug,...
                        'p',obj.p_val(2:end));
                    res_out = tnlp_results;
                    if tnlp_solver.stats.success
                        converged = true;
                    else
                        n_biactive = n_biactive + 1;
                    end 
                    %converged = true;
                    switch tnlp_solver.stats.return_status
                      case {'Solve_Succeeded', 'Solved_To_Acceptable_Level', 'Search_Direction_Becomes_Too_Small'}
                        converged = true;
                      otherwise
                        converged = false;
                    end
                end
                w_star = full(tnlp_results.x);
                w_star_orig = w_star;
                max(abs(w_init - w_star))
                [~,sort_idx] = sort(full(abs(w_init - w_star))); % sorted by difference from original stationary point

                nu = -full(tnlp_results.lam_g(ind_G));
                xi = -full(tnlp_results.lam_g(ind_H));
                
                % Index sets
                I_G = find(ind_0p+ind_00);
                I_H = find(ind_p0+ind_00);

                G_new = full(mpcc.G_fun(w_star_orig, obj.p_val(2:end)));
                H_new = full(mpcc.H_fun(w_star_orig, obj.p_val(2:end)));

                % multipliers (nu, xi already calculated above)
                lam_x_star = full(tnlp_results.lam_x);
                lam_g_star = full(tnlp_results.lam_g);
                %lam_x_star(ind_G) = 0; % Capturing multipliers for non-complementarity 
                %lam_x_star(ind_H) = 0; % box constraints
                lam_g_star = tnlp_results.lam_g;
                lam_g_star(ind_G) = 0; % Capturing multipliers for non-complementarity 
                lam_g_star(ind_H) = 0; % box constraints
                type_tol = a_tol^2;

                ind_00_new = G_new<a_tol & H_new<a_tol;

            end
            
            switch tnlp_solver.stats.return_status
              case {'Solve_Succeeded', 'Solved_To_Acceptable_Level', 'Search_Direction_Becomes_Too_Small'}
                
                if n_biactive
                    nu_biactive = nu(ind_00);
                    xi_biactive = xi(ind_00);
                    bound = 1.1*(max(abs([nu_biactive;xi_biactive]))+1e-10);

                    figure()
                    scatter(nu_biactive, xi_biactive, 50, 'o', 'LineWidth', 2);
                    xline(0,'k-.');
                    yline(0,'k-.');
                    xlim([-bound,bound]);
                    ylim([-bound,bound]);
                    
                    grid on;

                    if all(nu_biactive > -type_tol & xi_biactive > -type_tol)
                        stat_type = "S";
                        disp("Converged to S-stationary point")
                    elseif all((nu_biactive > -type_tol & xi_biactive > -type_tol) | (abs(nu_biactive.*xi_biactive) < type_tol))
                        stat_type = "M";
                        disp("Converged to M-stationary point")
                    elseif all(nu_biactive.*xi_biactive > -type_tol)
                        stat_type = "C";
                        disp("Converged to C-stationary point")
                    elseif all(nu_biactive > -type_tol | xi_biactive > -type_tol)
                        stat_type = "A";
                        disp("Converged to A-stationary point")
                    else
                        stat_type = "W";
                        disp("Converged to W-stationary point, or something has gone wrong")
                    end
                end
              otherwise
                stat_type = "?";
                disp("Could not converge to point from the end of homotopy");
            end
            
            % output tnlp results
            res_out = tnlp_results;
            polished_w = zeros(size(nlp.w));
            polished_w(nlp.ind_map) = w_star_orig;
            res_out.x = polished_w;
        end

        function [solution, improved_point, b_stat] = check_b_stationarity(obj, w, s_elastic)
            import casadi.*
            mpcc = obj.mpcc;
            nlp = obj.nlp;
            x0 = w(nlp.ind_map);
            f = mpcc.augmented_objective;
            x = mpcc.w; lbx = mpcc.lbw; ubx = mpcc.ubw;
            ind_g = setdiff(1:length(nlp.g),nlp.ind_g_comp);
            g = nlp.g(ind_g); lbg = nlp.lbg(ind_g); ubg = nlp.ubg(ind_g);
            G = mpcc.G_fun(mpcc.w, mpcc.p);
            H = mpcc.H_fun(mpcc.w, mpcc.p);
            G_old = full(mpcc.G_fun(x0, obj.p_val(2:end)));
            H_old = full(mpcc.H_fun(x0, obj.p_val(2:end)));

            if obj.solver_options.elasticity_mode == ElasticityMode.NONE
                a_tol = 1*sqrt(obj.solver_options.comp_tol);
            else
                a_tol = 1*sqrt(s_elastic);
            end
            ind_00 = G_old<a_tol & H_old<a_tol;
            n_biactive = sum(ind_00)
            ind_0p = G_old<H_old & ~ind_00;
            ind_p0 = H_old<G_old & ~ind_00;
            y_lpcc = H_old<G_old;
            
            lift_G = SX.sym('lift_G', length(G));
            lift_H = SX.sym('lift_H', length(H));
            g_lift = vertcat(lift_G - G, lift_H - H);

            x = vertcat(x, lift_G, lift_H);
            lblift_G = zeros(size(lift_G));
            ublift_G = inf*ones(size(lift_G));
            lblift_H = zeros(size(lift_G));
            ublift_H = inf*ones(size(lift_G));
            ind_mpcc = 1:length(lbx);
            ind_G = (1:length(lift_G))+length(lbx);
            ind_H = (1:length(lift_H))+length(lbx)+length(lift_G);

            lbx = vertcat(lbx, lblift_G, lblift_H);
            ubx = vertcat(ubx, ublift_G, ublift_H);
            G_init = G_old;
            %G_init(ind_00 | ind_0p) = 0;
            H_init = H_old;
            %H_init(ind_00 | ind_p0) = 0;
            x0 = vertcat(x0, G_init, H_init);

            g = vertcat(g,g_lift);
            lbg = vertcat(lbg,zeros(size(g_lift)));
            ubg = vertcat(ubg,zeros(size(g_lift)));
            
            p = mpcc.p;
            p0 = obj.p_val(2:end);

            solver_settings = struct();
            solver_settings.max_iter = 20;
            solver_settings.tol = 1e-6;
            solver_settings.Delta_TR_init = 10;
            solver_settings.Delta_TR_min = 1e-4;
            solver_settings.verbose_solver = 1;
            solver_settings.tighten_bounds_in_lpcc = false;
            solver_settings.BigM = 1e2;

            if ~n_biactive
                solver_settings.fixed_y_lpcc = y_lpcc;
            end
            %% The problem
            nlp = struct('x', x, 'f', f, 'g', g, 'comp1', lift_G, 'comp2', lift_H, 'p', p);
            solver_initalization = struct('x0', x0, 'lbx', lbx, 'ubx', ubx,'lbg', lbg, 'ubg', ubg, 'p', p0, 'y_lpcc', y_lpcc);
            solution = bStationarityOracle(nlp,solver_initalization,solver_settings);
            cpu_time_filterSMPCC2 = toc;

            improved_point = solution.x(ind_mpcc);

            b_stat = ~solution.oracle_status;
        end
    end

    methods(Access=private)
        function obj = construct_relaxation_homotopy(obj)
            import casadi.*
            mpcc = obj.mpcc;
            solver_options = obj.solver_options;
            tic;
            % calculate homotopy
            if solver_options.N_homotopy == 0
                solver_options.N_homotopy = ceil(abs(log(solver_options.sigma_N / solver_options.sigma_0) / log(solver_options.homotopy_update_slope)));
                % TODO: compute
                if ~strcmp(solver_options.homotopy_update_rule, 'linear')
                    warning('computing N_homotopy automatically only supported for linear homotopy_update_rule');
                end
            end

            if isstruct(obj.mpcc)
                nlp = obj.generate_generic_nlp();
            else
                nlp = NosnocNLP(solver_options, mpcc);
            end

            obj.nlp = nlp;
            w = nlp.w;
            g = nlp.g;
            p = nlp.p;
            if ~isempty(solver_options.ipopt_callback)
                solver_options.opts_casadi_nlp.iteration_callback = NosnocIpoptCallback('a_callback', model, nlp, solver_options, length(w),length(g),length(p));
            end
            
            obj.solver = obj.plugin.construct_solver(nlp, solver_options);

            if ~isempty(solver_options.ipopt_callback)
                solver_options.opts_casadi_nlp.iteration_callback.solver = solver;
            end

            if solver_options.print_level > 5
                nlp.print();
            end

            solver_generating_time = toc;
            if solver_options.print_level >=2
                fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
            end
        end

        function obj = construct_direct(obj)
            error('TODO: not yet implemented')
        end

        function [results, stats] = solve_relaxation_homotopy(obj)
            import casadi.*;
            solver = obj.solver;
            mpcc = obj.mpcc;
            solver_options = obj.solver_options;
            nlp = obj.nlp;
            plugin = obj.plugin;

            comp_res = mpcc.comp_res;

            % Initial conditions
            sigma_k = solver_options.sigma_0;
            w0 = nlp.w0;
            w0_mpcc = w0;
            w0_mpcc(nlp.ind_elastic) = [];
            lbw = nlp.lbw; ubw = nlp.ubw;
            lbg = nlp.lbg; ubg = nlp.ubg;
            obj.compute_initial_parameters();

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.wall_time = [];
            stats.cpu_time_total = 0;
            stats.wall_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res(w0_mpcc, obj.p_val(2:end)))];


            % Initialize Results struct
            results = struct;
            results.W = w0;
            results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;
            last_iter_failed = 0;
            timeout = 0;

            if solver_options.print_level >= 3
                plugin.print_nlp_iter_header();
            end

            while (abs(complementarity_iter) > solver_options.comp_tol || last_iter_failed) &&...
                    ii < solver_options.N_homotopy &&...
                    (sigma_k > solver_options.sigma_N || ii == 0) &&...
                    ~timeout
                % homotopy parameter update
                last_iter_failed = 0;
                if ii == 0
                    sigma_k = solver_options.sigma_0;
                else
                    if isequal(solver_options.homotopy_update_rule,'linear')
                        sigma_k = solver_options.homotopy_update_slope*sigma_k;
                    elseif isequal(solver_options.homotopy_update_rule,'superlinear')
                        sigma_k = max(solver_options.sigma_N,min(solver_options.homotopy_update_slope*sigma_k,sigma_k^solver_options.homotopy_update_exponent));
                    else
                        error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
                    end
                end
                stats.sigma_k = [stats.sigma_k, sigma_k];
                obj.p_val(1) = sigma_k;

                % Using multi-solver
                if iscell(solver)
                    start = cputime;
                    tic;
                    nlp_results = solver{ii+1}('x0', w0,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'p',obj.p_val);
                    cpu_time_iter = cputime - start;
                    wall_time_iter = toc;
                    stats.solver_stats = [stats.solver_stats, solver{ii+1}.stats];
                else
                    if ii ~= 0

                        % TODO(Anton) Lets push on casadi devs to allow for changing of options after construction
                        %             (or maybe do it ourselves) and then remove this hack
                        if obj.solver_options.timeout_cpu
                            solver = plugin.construct_solver(nlp, solver_options, solver_options.timeout_cpu - stats.cpu_time_total);
                        elseif obj.solver_options.timeout_wall
                            solver = plugin.construct_solver(nlp, solver_options, solver_options.timeout_wall - stats.wall_time_total);
                        end % HACK ENDS HERE
                        
                        start = cputime;
                        tic;
                        nlp_results = solver('x0', w0,...
                            'lbx', lbw,...
                            'ubx', ubw,...
                            'lbg', lbg,...
                            'ubg', ubg,...
                            'lam_x0', nlp_results.lam_x,...
                            'lam_g0', nlp_results.lam_g,...
                            'p',obj.p_val);
                            %'lam_x0', nlp_results.lam_x,...
                            %'lam_g0', nlp_results.lam_g);
                        cpu_time_iter = cputime - start;
                        wall_time_iter = toc;
                    else
                        start = cputime;
                        tic;
                        nlp_results = solver('x0', w0,...
                            'lbx', lbw,...
                            'ubx', ubw,...
                            'lbg', lbg,...
                            'ubg', ubg,...
                            'p',obj.p_val);
                        cpu_time_iter = cputime - start;
                        wall_time_iter = toc;
                    end
                    solver_stats = plugin.cleanup_solver_stats(solver.stats);
                    
                    stats.solver_stats = [stats.solver_stats, solver_stats];
                end

                results.nlp_results = [results.nlp_results, nlp_results];

                last_iter_failed = plugin.check_iteration_failed(stats);
                timeout = plugin.check_timeout(stats);

                if last_iter_failed
                    obj.print_infeasibility();
                end

                % update timing stats
                stats.cpu_time = [stats.cpu_time,cpu_time_iter];
                stats.cpu_time_total = stats.cpu_time_total + cpu_time_iter;
                stats.wall_time = [stats.wall_time, wall_time_iter];
                stats.wall_time_total = stats.wall_time_total + wall_time_iter;

                if solver_options.timeout_cpu && (stats.cpu_time_total > solver_options.timeout_cpu)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                if solver_options.timeout_wall && (stats.wall_time_total > solver_options.timeout_wall)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                % update results output.
                w_opt = plugin.w_opt_from_results(nlp_results);
                w_opt_mpcc = w_opt;
                w_opt_mpcc(nlp.ind_elastic) = [];
                results.W = [results.W,w_opt]; % all homotopy iterations
                w0 = w_opt;

                % update complementarity and objective stats
                complementarity_iter = full(comp_res(w_opt_mpcc, obj.p_val(2:end)));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(obj.nlp.objective_fun(w_opt, obj.p_val));
                stats.objective = [stats.objective, objective];

                % update counter
                ii = ii+1;
                % Verbose
                if solver_options.print_level >= 3
                    plugin.print_nlp_iter_info(stats)
                end
            end

            % polish homotopy solution with fixed active set.
            % TODO fix this!
            if solver_options.polishing_step
                [results] = polish_homotopy_solution(model,problem,solver_options,results,sigma_k);
                complementarity_iter = results.complementarity_iter;
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                W = [W,results.w];
            end

            if solver_options.calculate_stationarity_type
                if last_iter_failed || ~obj.complementarity_tol_met(stats) || timeout
                    stat_type = "?";
                    disp("Not checking stationarity due to failure of homotopy to converge");
                else
                    [polished_w, res_out, stat_type, n_biactive] = obj.calculate_stationarity(results.nlp_results(end), false);
                end
                if stat_type ~= "?"
                    results.W = [results.W,polished_w];
                    results.nlp_results = [results.nlp_results, res_out];
                    complementarity_iter = full(comp_res(polished_w(nlp.ind_map), obj.p_val(2:end)));
                    stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                end
                stats.stat_type = stat_type;
            end
            
            % number of iterations
            stats.homotopy_iterations = ii;

            results = obj.extract_results_nlp(results);

            % check if solved to required accuracy
            stats.converged = obj.complementarity_tol_met(stats) && ~last_iter_failed && ~timeout;
            stats.constraint_violation = obj.compute_constraint_violation(results.w, nlp.g_fun);
        end

        function obj = solve_direct(obj)
            error('TODO: not yet implemented')
        end

        function results = extract_results_nlp(obj, results)
            import casadi.*
            mpcc = obj.mpcc;
            plugin = obj.plugin;
            % Store differential states
            w_opt = plugin.w_opt_from_results(results.nlp_results(end));
            results.w = w_opt;

            w_mpcc = w_opt;
            w_mpcc(obj.nlp.ind_elastic) = [];

            names = get_result_names_from_settings(obj.mpcc.problem_options);
            % populate outputs
            for name=names
                results = form_structured_output(mpcc, w_mpcc, name, results);
            end

            % handle x0 properly
            x0 = w_mpcc(mpcc.ind_x0);
            results.x = [x0, results.x];
            results.extended.x = [x0, results.extended.x];

            if ~isempty(mpcc.ind_v_global)
                v_global = w_mpcc(mpcc.ind_v_global);
                results.v_global = v_global;
            end
            
            u = w_mpcc([mpcc.ind_u{:}]);
            u = reshape(u,obj.mpcc.model.dims.n_u,obj.mpcc.problem_options.N_stages);

            results.u = u;

            if mpcc.problem_options.time_optimal_problem
                T_opt = w_mpcc(mpcc.ind_t_final);
            else
                T_opt = mpcc.problem_options.T;
            end
            results.T = T_opt;

            if obj.mpcc.problem_options.use_fesd
                h_opt = w_mpcc(flatten_ind(mpcc.ind_h))';
            else
                h_opt = [];
                % if settings.time_optimal_problem && ~settings.use_speed_of_time_variables
                %     T = T_opt;
                % end
                for ii = 1:obj.mpcc.problem_options.N_stages
                    h_opt = [h_opt,obj.mpcc.problem_options.T/(obj.mpcc.problem_options.N_stages*obj.mpcc.problem_options.N_finite_elements(ii))*ones(1, obj.mpcc.problem_options.N_finite_elements(ii))];
                end
            end
            results.h = h_opt;

            t_grid = cumsum([0,h_opt]);

            if obj.mpcc.problem_options.use_speed_of_time_variables
                s_sot = w_mpcc(flatten_ind(mpcc.ind_sot));
                if ~obj.mpcc.problem_options.local_speed_of_time_variable
                    s_sot = s_sot*ones(obj.mpcc.problem_options.N_stages,1);
                end
                results.s_sot = s_sot;
            end
            
            %% Adapt the grid in case of time optimal problems
            if obj.mpcc.problem_options.time_optimal_problem
                if obj.mpcc.problem_options.use_speed_of_time_variables
                    s_sot = w_mpcc(flatten_ind(mpcc.ind_sot));
                    if ~obj.mpcc.problem_options.local_speed_of_time_variable
                        s_sot = s_sot*ones(obj.mpcc.problem_options.N_stages,1);
                    end
                    results.s_sot = s_sot;
                    h_rescaled = [];
                    ind_prev = 1;
                    for ii = 1:obj.mpcc.problem_options.N_stages
                        h_rescaled = [h_rescaled,h_opt(ind_prev:obj.mpcc.problem_options.N_finite_elements(ii)+ind_prev-1).*s_sot(ii)];
                        ind_prev = ind_prev+obj.mpcc.problem_options.N_finite_elements(ii);
                    end
                    t_grid = cumsum([0,h_rescaled]);
                else
                    t_grid = cumsum([0,h_opt]);
                end
            end
            ind_t_grid_u = cumsum([1; obj.mpcc.problem_options.N_finite_elements]);

            if obj.mpcc.problem_options.dcs_mode == DcsMode.CLS
                x_with_impulse = x0;
                t_with_impulse = kron(t_grid, ones(2,1));
                for ii=1:size(results.structured.x,1)
                    for jj=1:size(results.structured.x,2)
                        x_with_impulse = [x_with_impulse,results.structured.x_left_bp{ii,jj}];
                        x_fe = results.structured.x{ii,jj};
                        x_with_impulse = [x_with_impulse, x_fe(:,end)];
                    end
                end
                results.x_with_impulse = x_with_impulse;
                results.t_with_impulse = t_with_impulse(1:end-1);
            end


            results.t_grid = t_grid;
            results.t_grid_u = t_grid(ind_t_grid_u);

            results.f = full(results.nlp_results(end).f);
            results.g = full(results.nlp_results(end).g);
            results.objective = full(obj.nlp.objective_fun(w_opt,obj.p_val));
        end

        function nlp = generate_generic_nlp(obj)
            import casadi.*
            solver_options = obj.solver_options;
            mpcc = obj.mpcc;

            % This is an ugly hack and the backend rewrite should perhaps handle this smarter
            nlp = struct();
            w = mpcc.w;
            lbw = mpcc.lbw;
            ubw = mpcc.ubw;
            w0 = mpcc.w0;
            g = mpcc.g;
            lbg = mpcc.ubg;
            ubg = mpcc.lbg;
            p = mpcc.p;
            p0 = mpcc.p0;
            G = mpcc.G;
            H = mpcc.H;
            augmented_objective = mpcc.f;

            casadi_symbolic_mode = ['casadi.' w.type_name()];

            psi_fun = obj.solver_options.psi_fun;
            sigma_p = define_casadi_symbolic(casadi_symbolic_mode, 'sigma_p', 1);
            p = [sigma_p;p];
            p0 = [1;p0]
            if obj.solver_options.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            elseif obj.solver_options.elasticity_mode == ElasticityMode.ELL_INF
                sigma = define_casadi_symbolic(casadi_symbolic_mode, 's_elastic', 1);
                w = [w;sigma];
                lbw = [lbw;obj.solver_options.s_elastic_min];
                ubw = [ubw;obj.solver_options.s_elastic_max];
                w0 = [w0;obj.solver_options.s_elastic_0];
                if obj.solver_options.elastic_scholtes
                    obj.solver_options.s_elastic_max = inf;
                    g = [g;sigma - sigma_p];
                    lbg = [lbg;-inf];
                    ubg = [ubg;0];
                end
                if solver_options.objective_scaling_direct
                    augmented_objective = augmented_objective + (1/sigma_p)*sigma;
                else
                    augmented_objective = sigma_p*augmented_objective + sigma;
                end
            else
                n_comp = size(mpcc.G, 1);
                sigma = define_casadi_symbolic(casadi_symbolic_mode, 's_elastic', n_comp);
                w = [w,sigma];
                lbw = [lbw;obj.solver_options.s_elastic_min*ones(n_comp,1)];
                ubw = [ubw;obj.solver_options.s_elastic_max*ones(n_comp,1)];
                w0 = [w0;obj.solver_options.s_elastic_0*ones(n_comp,1)];
                if obj.solver_options.elastic_scholtes
                    obj.solver_options.s_elastic_max = inf;
                    g = [g;sigma - sigma_p];
                    lbg = [lbg;-inf*ones(n_comp,1)];
                    ubg = [ubg;0*ones(n_comp,1)];
                end
                if solver_options.objective_scaling_direct
                    augmented_objective = augmented_objective + (1/sigma_p)*sum1(sigma);
                else
                    augmented_objective = sigma_p*augmented_objective + sum1(sigma);
                end
            end

            expr = psi_fun(mpcc.G, mpcc.H, sigma);
            [lb, ub, expr] = generate_mpcc_relaxation_bounds(expr, obj.solver_options);
            g = [g;expr];
            lbg = [lbg;lb];
            ubg = [ubg;ub];

            nlp.w = w;
            nlp.lbw = lbw;
            nlp.ubw = ubw;
            nlp.w0 = w0;
            nlp.g = g;
            nlp.lbg = lbg;
            nlp.ubg = ubg;
            nlp.p = p;
            nlp.p0 = p0;
            nlp.augmented_objective = augmented_objective;
        end

        function [results,stats] = solve_generic_homotopy(obj)
            import casadi.*;
            solver = obj.solver;
            mpcc = obj.mpcc;
            solver_options = obj.solver_options;
            nlp = obj.nlp;
            plugin = obj.plugin;

            comp_res = Function('comp_res', {nlp.w,nlp.p}, {mmax(mpcc.G.*mpcc.H)});
            objective_fun = Function('objective', {mpcc.w, nlp.p}, {mpcc.f});

            % Initial conditions
            sigma_k = solver_options.sigma_0;
            w0 = nlp.w0;
            lbw = nlp.lbw; ubw = nlp.ubw;
            lbg = nlp.lbg; ubg = nlp.ubg;

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.wall_time = [];
            stats.cpu_time_total = 0;
            stats.wall_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res(w0, nlp.p0))];


            % Initialize Results struct
            results = struct;
            results.W = w0;
            results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;
            last_iter_failed = 0;
            timeout = 0;

            if solver_options.print_level >= 3
                plugin.print_nlp_iter_header();
            end

            while (abs(complementarity_iter) > solver_options.comp_tol || last_iter_failed) &&...
                    ii < solver_options.N_homotopy &&...
                    (sigma_k > solver_options.sigma_N || ii == 0) &&...
                    ~timeout
                % homotopy parameter update
                last_iter_failed = 0;
                if ii == 0
                    sigma_k = solver_options.sigma_0;
                else
                    if isequal(solver_options.homotopy_update_rule,'linear')
                        sigma_k = solver_options.homotopy_update_slope*sigma_k;
                    elseif isequal(solver_options.homotopy_update_rule,'superlinear')
                        sigma_k = max(solver_options.sigma_N,min(solver_options.homotopy_update_slope*sigma_k,sigma_k^solver_options.homotopy_update_exponent));
                    else
                        error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
                    end
                end
                stats.sigma_k = [stats.sigma_k, sigma_k];
                obj.p_val(1) = sigma_k;

                if ii ~= 0

                    % TODO(Anton) Lets push on casadi devs to allow for changing of options after construction
                    %             (or maybe do it ourselves) and then remove this hack
                    if obj.solver_options.timeout_cpu
                        solver = plugin.construct_solver(nlp, solver_options, solver_options.timeout_cpu - stats.cpu_time_total);
                    elseif obj.solver_options.timeout_wall
                        solver = plugin.construct_solver(nlp, solver_options, solver_options.timeout_wall - stats.wall_time_total);
                    end % HACK ENDS HERE
                    
                    start = cputime;
                    tic;
                    nlp_results = solver('x0', w0,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'lam_x0', nlp_results.lam_x,...
                        'lam_g0', nlp_results.lam_g,...
                        'p',obj.p_val);
                    %'lam_x0', nlp_results.lam_x,...
                    %'lam_g0', nlp_results.lam_g);
                    cpu_time_iter = cputime - start;
                    wall_time_iter = toc;
                else
                    start = cputime;
                    tic;
                    nlp_results = solver('x0', w0,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'p',obj.p_val);
                    cpu_time_iter = cputime - start;
                    wall_time_iter = toc;
                end
                solver_stats = plugin.cleanup_solver_stats(solver.stats);
                
                stats.solver_stats = [stats.solver_stats, solver_stats];

                results.nlp_results = [results.nlp_results, nlp_results];

                last_iter_failed = plugin.check_iteration_failed(stats);
                timeout = plugin.check_timeout(stats);

                if last_iter_failed
                    obj.print_infeasibility();
                end

                % update timing stats
                stats.cpu_time = [stats.cpu_time,cpu_time_iter];
                stats.cpu_time_total = stats.cpu_time_total + cpu_time_iter;
                stats.wall_time = [stats.wall_time, wall_time_iter];
                stats.wall_time_total = stats.wall_time_total + wall_time_iter;

                if solver_options.timeout_cpu && (stats.cpu_time_total > solver_options.timeout_cpu)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                if solver_options.timeout_wall && (stats.wall_time_total > solver_options.timeout_wall)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                % update results output.
                w_opt = plugin.w_opt_from_results(nlp_results);
                w_opt_mpcc = w_opt(1:length(mpcc.w));
                results.W = [results.W,w_opt]; % all homotopy iterations
                w0 = w_opt;

                % update complementarity and objective stats
                complementarity_iter = full(comp_res(w_opt, obj.p_val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(objective_fun(w_opt_mpcc, obj.p_val));
                stats.objective = [stats.objective, objective];

                % update counter
                ii = ii+1;
                % Verbose
                if solver_options.print_level >= 3
                    plugin.print_nlp_iter_info(stats)
                end
            end

            results.w = w_opt;
            
            % number of iterations
            stats.homotopy_iterations = ii;

            % check if solved to required accuracy
            stats.converged = obj.complementarity_tol_met(stats) && ~last_iter_failed && ~timeout;
            g_fun = Function('g_fun', {nlp.w, nlp.p}, {nlp.g});
            stats.constraint_violation = obj.compute_constraint_violation(results.w, g_fun);
        end
    end
end
