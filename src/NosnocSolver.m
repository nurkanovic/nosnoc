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
        solver % Nosnoc solver function
    end

    methods
        function obj = NosnocSolver(mpcc, solver_options)
            import casadi.*
            obj.mpcc = mpcc;
            obj.solver_options = solver_options;
            opts = struct;
            opts.relaxation = solver_options;
            mpcc_struct = struct;
            mpcc_struct.x = mpcc.w;
            mpcc_struct.g = mpcc.g;
            mpcc_struct.p = mpcc.p;
            mpcc_struct.f = mpcc.augmented_objective;
            mpcc_struct.G = mpcc.cross_comps(:,1);
            mpcc_struct.H = mpcc.cross_comps(:,2);
            
            obj.solver = nosnoc.solver.mpccsol('nosnoc_solver', 'relaxation', mpcc_struct, opts);
        end

        function set(obj, type, val)
            mpcc = obj.mpcc;
            if strcmp(type, 'w0')
                if length(obj.mpcc.w0) == length(val)
                    obj.mpcc.w0 = val;
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
                                    obj.mpcc.w0(v{1}) = val{ii};
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
                                        obj.mpcc.w0(v{1}) = val{ii,jj};
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
                        obj.mpcc.w0(flat_ind) = val;
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end

        function [results,stats] = solve(obj)
            mpcc = obj.mpcc;
            p0 = [mpcc.p0;mpcc.compute_initial_parameters(mpcc.w0(mpcc.ind_x0))]; % TODO(@anton) correct?
            results = obj.solver('x0', mpcc.w0,...
                'lbx', mpcc.lbw,...
                'ubx', mpcc.ubw,...
                'lbg', mpcc.lbg,...
                'ubg', mpcc.ubg,...
                'p', p0);
            stats = obj.solver.stats;
            results = obj.extract_results(results);
        end
        % TODO(@anton) fix these two or more likely move them to nosnoc.solver.RelaxationSolver
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

        % TODO(@anton) most of this will die and be replaced with vdx :)
        function results = extract_results(obj, results)
            import casadi.*
            mpcc = obj.mpcc;
            % Store differential states
            w_opt = full(results.x);
            results.w = w_opt;
            w_mpcc = w_opt;

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

            results.f = full(results.f);
            results.g = full(results.g);
            results.objective = results.f;
        end
    end
end
