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

classdef MpccSolver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc % Either a struct with the (possibly optional) fields (f, p, w, g, G, H) or a subclass of vdx.problems.Mpcc.
        nlp % The relaxed/smoothed nlp which is solved in a homotopy loop with a decreasing relaxation/smoothing parameter. 
        opts % Options object.
        stats % Struct with stats of the last solve.
        nlp_solver
        plugin % NLP solver: Ipopt, Snopt, Worhop or Uno.
        relaxation_type
    end

    properties (Access=private)
        ind_scalar_G
        ind_scalar_H
        ind_nonscalar_G
        ind_nonscalar_H
        ind_map_G % Map from lifted indices to original indices in G.
        ind_map_H % Map from lifted indices to original indices in H.

        G_fun
        H_fun
        comp_res_fun
        f_mpcc_fun
    end

    methods (Access=public)
        
        function obj=MpccSolver(relaxation_type, mpcc, opts)
            import casadi.*
            import nosnoc.solver.*
            casadi_symbolic_mode = split(class(mpcc.x), '.');
            casadi_symbolic_mode = casadi_symbolic_mode{end};
            % preprocess empty fields
            if ~isfield(mpcc, 'g') || isempty(mpcc.g)
                mpcc.g = casadi.(casadi_symbolic_mode)([]); 
            end
            if ~isfield(mpcc, 'p') || isempty(mpcc.p)
                mpcc.p = casadi.(casadi_symbolic_mode)([]); 
            end
            obj.relaxation_type = relaxation_type;
            obj.mpcc = mpcc;
            obj.opts = opts;
            opts.preprocess();
            stats = struct;
            
            if isa(mpcc, 'vdx.problems.Mpcc') % We can properly interleave complementarities if we get a vdx.Mpcc
                use_vdx = true;
                % TODO(@anton) implement this
            else % Otherwise use vdx internally anyway but be sad about interleaving
                use_vdx = false;

                % create an nlp (in the form of a vdx.Problem) from the mpcc data (f,x,p,g).
                % We track the indices of this data in vdx.Variables: mpcc_w, mpcc_p, mpcc_g.
                nlp = vdx.Problem('casadi_type', casadi_symbolic_mode);
                nlp.w.mpcc_w = {mpcc.x}; % nlp.w.mpcc_w is a vdx.Variable (depth 0) which stores mpcc optimization variables
                nlp.p.mpcc_p = {mpcc.p}; % nlp.p.mpcc_p is a vdx.Variable (depth 0) which stores mpcc parameters
                nlp.g.mpcc_g = {mpcc.g}; % nlp.g.mpcc_g is a vdx.Variable (depth 0) which stores mpcc general constraints
                nlp.p.sigma_p = {{'sigma_p', 1}, opts.sigma_0}; % sigma_p is the homotopy parameter (a scalar), with the value sigma_0 from the options.
                nlp.f = mpcc.f;
                n_c = length(mpcc.G); % number of complementarity constraints

                % Create relaxation slacks/parameters
                switch opts.homotopy_steering_strategy
                  case HomotopySteeringStrategy.DIRECT
                    % nlp.p.sigma_p(): sigma is a parameter/variable that has no indices
                    sigma = nlp.p.sigma_p(); 
                  case HomotopySteeringStrategy.ELL_INF
                    % adding a scalar elastic variable to nlp.w which augments the original mpcc.w
                    nlp.w.s_elastic = {{'s_elastic', 1}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                    sigma = nlp.w.s_elastic(); % Here s_elastic takes the role of sigma in direct, and sigma_p is used to define a penalty parameter for the elastic variable s_elastic
                    if opts.objective_scaling_direct
                        nlp.f = nlp.f + (1/nlp.p.sigma_p())*sigma; % penalize the elastic more and more with decreasing sigma_p
                    else
                        nlp.f = nlp.p.sigma_p()*nlp.f + sigma; % reduce the weight of the initial objective with decreasing sigma_p
                    end
                  case HomotopySteeringStrategy.ELL_1
                    % adding elastic variables to nlp.w which augments the original mpcc.w
                    % Remark: ELL_1 with s_elastic is equivalent to usually Ell_1 penality approach, but this indirect way helps
                    % to add some constraint on s_elastic (which  avoids unbounded problems somtimes, and it can also improve convergence)
                    nlp.w.s_elastic = {{'s_elastic', n_c}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                    sigma = nlp.w.s_elastic();
                    sum_elastic = sum1(sigma);
                    if opts.objective_scaling_direct
                        nlp.f = nlp.f + (1/nlp.p.sigma_p())*sum_elastic;
                    else
                        nlp.f = nlp.p.sigma_p()*nlp.f + sum_elastic;
                    end
                end

                [ind_scalar_G,ind_nonscalar_G, ind_map_G] = find_nonscalar(mpcc.G,mpcc.x);
                [ind_scalar_H,ind_nonscalar_H, ind_map_H] = find_nonscalar(mpcc.H,mpcc.x);
                obj.ind_nonscalar_G = ind_nonscalar_G;
                obj.ind_nonscalar_H = ind_nonscalar_H;
                obj.ind_scalar_G = ind_scalar_G;
                obj.ind_scalar_H = ind_scalar_H;
                obj.ind_map_G = ind_map_G;
                obj.ind_map_H = ind_map_H;
                % possibly lift complementarities
                if opts.lift_complementarities
                    nlp.w.G_lift = {{'G', length(ind_nonscalar_G)}, 0, inf};
                    G = casadi.(casadi_symbolic_mode)(size(mpcc.H,1), 1);
                    G(ind_scalar_G) = mpcc.G(ind_scalar_G);
                    G(ind_nonscalar_G) = nlp.w.G_lift();

                    nlp.w.H_lift = {{'H', length(ind_nonscalar_H)}, 0, inf};
                    H = casadi.(casadi_symbolic_mode)(size(mpcc.H,1), 1);
                    H(ind_scalar_H) = mpcc.H(ind_scalar_H);
                    H(ind_nonscalar_H) = nlp.w.H_lift();
                    
                    nlp.g.G_lift = {mpcc.G(ind_nonscalar_G)-G(ind_nonscalar_G)};
                    nlp.g.H_lift = {mpcc.H(ind_nonscalar_H)-H(ind_nonscalar_H)};    
                else
                    G = mpcc.G;
                    H = mpcc.H;
                end

                % apply relaxation (defines a particular method, e.g. Scholtes, Kanzow-Schwartz, Fischer Burmeister, etc.)
                psi_fun = get_psi_fun(MpccMethod(obj.relaxation_type), opts.normalize_homotopy_update);
                lb = [];
                ub = [];
                g_comp_expr = [];
                n_comp_pairs = size(G, 1);
                for ii=1:n_comp_pairs
                    g_comp_expr_i = psi_fun(G(ii), H(ii), sigma);
                    [lb_i, ub_i, g_comp_expr_i] = generate_mpcc_relaxation_bounds(g_comp_expr_i, obj.relaxation_type);
                    lb = [lb;lb_i];
                    ub = [ub;ub_i];
                    g_comp_expr = [g_comp_expr;g_comp_expr_i];
                end
                nlp.g.complementarities(0) = {g_comp_expr, lb, ub};

                if ~opts.assume_lower_bounds && ~opts.lift_complementarities % Lower bounds on G, H, not already present in MPCC
                    if ~isempty(obj.ind_nonscalar_G)
                        nlp.g.G_lower_bounds = {mpcc.G(obj.ind_nonscalar_G), 0, inf};
                    end
                    if ~isempty(obj.ind_nonscalar_H)
                        nlp.g.H_lower_bounds = {mpcc.H(obj.ind_nonscalar_H), 0, inf};
                    end
                end
                % Get nlpsol plugin
                switch opts.solver
                  case 'ipopt'
                    obj.plugin = nosnoc.solver.plugins.Ipopt();
                  case 'snopt'
                    obj.plugin = nosnoc.solver.plugins.Snopt();
                  case 'worhp'
                    obj.plugin = nosnoc.solver.plugins.Worhp();
                  case 'uno'
                    obj.plugin = nosnoc.solver.plugins.Uno();
                end

                % Construct solver
                obj.plugin.construct_solver(nlp, opts);
                obj.nlp = nlp;
            end
        end

        function met = complementarity_tol_met(obj, stats)
            last_stats = stats.solver_stats(end);
            met = 0;
            if abs(stats.complementarity_stats(end)) < obj.opts.complementarity_tol + eps(obj.opts.complementarity_tol) % add epsilon to check.
                met = 1;
            end
        end

        function violation = compute_constraint_violation(obj, w, g)
            nlp = obj.nlp;
            % We get the maximum violation of g: (max(max(lbg - g(w,p), 0), max(g(w,p) - ubg, 0))),
            % and w: (max(max(lbw - w, 0), max(w - ubw, 0)))
            % bounds from the last solve of the relaxed nlp.
            violation = max([nlp.g.mpcc_g().violation; nlp.w.mpcc_w().violation]);
        end

        function out = cat(dim,varargin)
            error('Concatenation not supported.')
        end

        function varargout = size(obj,varargin)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            varargout = {1};
        end
        
        function ind = end(obj,k,n)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            ind = 1;
        end

        function [w_polished, res_out, stat_type, n_biactive] = check_multiplier_based_stationarity(obj, x0, exitfast, complementarity_constraints_lifted)
        % exitfast: Exit without solving the TNLP if we have no biactive constraints (point must be S-stationary).
        % complementarity_constraints_lifted: Additionally lift necessary complementariy constraints when necessary.
            import casadi.*
            if ~exist('exitfast', 'var')
                exitfast = false;
            end
            if ~exist('complementarity_constraints_lifted', 'var')
                complementarity_constraints_lifted = true;
            end
            stat_type = "?";
            nlp = obj.nlp;
            mpcc = obj.mpcc;

            % This method generates a TNLP (tightened NLP) using the results of a homotopy solver.
            % In order to do so we first identify the active set of the complementarity constraints which involves identifying which
            % branch of the complementarity each element of G(w) and H(w) are, or whether the solution is biactive at that point.
            %
            % This algorithm also attempts to solve multiple TNLPs by iteratively adding ambiguous complementarity pairs to the bi-active set.
            % We do this starting from all the possibly biactive points 
            % and attempting to solve the corresponding TNLP. We then remove the ambiguous complementarities in order of inf-norm distance from the origin in G-H space.
            % This algorithm terminates when the correct biactive set is found (a TNLP converges to a point close enough to the solution of the homotopy), or until
            % all ambiguous complementarity pairs are exhausted, at which point we fail to calculate the stationarity type.
            %
            % For a more in depth explanation see the PhD theses of Armin Nurkanovic (Section 2.3, https://publications.syscop.de/Nurkanovic2023f.pdf)
            % or Alexandra Schwartz (Section 5.2, https://opus.bibliothek.uni-wuerzburg.de/opus4-wuerzburg/frontdoor/index/index/docId/4977).
            % A concise overview can also be found in https://arxiv.org/abs/2312.11022.
            
            w_orig = x0;
            lam_x = nlp.w.mpcc_w().mult;
            lam_g = nlp.g.mpcc_g().mult;

            G = mpcc.G;
            H = mpcc.H;
            G_old = full(obj.G_fun(w_orig, nlp.p.mpcc_p().val));
            H_old = full(obj.H_fun(w_orig, nlp.p.mpcc_p().val));

            % We take the square root of the complementarity tolerance for the activity tolerance as it is the upper bound for G, H, for any of the smoothing approaches.
            % i.e. in the worst case G=sqrt(complementarity_tol) and H=sqrt(comp_tol).
            % We add 1e-14 (a value near machine precision) in order to accurately identify points on the boundary.
            if obj.opts.homotopy_steering_strategy == HomotopySteeringStrategy.DIRECT
                a_tol = sqrt(obj.opts.complementarity_tol) + 1e-14;
            else
                a_tol = sqrt(max(nlp.w.s_elastic().res)) + 1e-14;
            end
            
            ind_00 = G_old<a_tol & H_old<a_tol;
            n_biactive = sum(ind_00);
            if exitfast && n_biactive == 0
                w_polished = [];
                res_out = [];
                stat_type = "S";
                disp("Converged to S-stationary point")
                return
            end
            min_dists = max(G_old, H_old);
            [~, min_idx] = sort(min_dists);

            if complementarity_constraints_lifted
                w = nlp.w.mpcc_w(); lbw_orig = nlp.w.mpcc_w().lb; ubw_orig = nlp.w.mpcc_w().ub;
                g = nlp.g.mpcc_g(); lbg = nlp.g.mpcc_g().lb; ubg = nlp.g.mpcc_g().ub;
                p = nlp.p.mpcc_p();
                f = mpcc.f;
                
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
                n_biactive = n_biactive;
                while ~converged && n_biactive >= 0
                    idx_00 = min_idx(1:n_biactive);
                    ind_00 = false(length(G),1);
                    ind_00(idx_00) = true;
                    
                    ind_0p = G_old<a_tol & ~ind_00;
                    ind_p0 = H_old<a_tol & ~ind_00;

                    lblift_G = -inf*ones(size(lift_G));
                    ublift_G = inf*ones(size(lift_G));
                    lblift_H = -inf*ones(size(lift_G));
                    ublift_H = inf*ones(size(lift_G));
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
                    opts_casadi_nlp.ipopt.tol = obj.opts.complementarity_tol;
                    opts_casadi_nlp.ipopt.acceptable_tol = sqrt(obj.opts.complementarity_tol);
                    opts_casadi_nlp.ipopt.acceptable_dual_inf_tol = sqrt(obj.opts.complementarity_tol);
                    opts_casadi_nlp.ipopt.fixed_variable_treatment = 'relax_bounds';
                    opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
                    opts_casadi_nlp.ipopt.bound_relax_factor = 1e-12;
                    opts_casadi_nlp.ipopt.linear_solver = obj.opts.opts_casadi_nlp.ipopt.linear_solver;
                    opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
                    opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
                    opts_casadi_nlp.ipopt.nlp_scaling_method = 'none';
                    default_tol = 1e-7;
                    opts_casadi_nlp.ipopt.tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.resto_failure_feasibility_threshold = 0;
                    opts_casadi_nlp.ipopt.print_level = obj.opts.opts_casadi_nlp.ipopt.print_level;
                    opts_casadi_nlp.print_time = obj.opts.opts_casadi_nlp.print_time;
                    opts_casadi_nlp.ipopt.sb = 'yes';
                    
                    
                    tnlp_solver = nlpsol('tnlp', 'ipopt', casadi_nlp, opts_casadi_nlp);
                    tnlp_results = tnlp_solver('x0', w_init,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'lam_g0', lam_g_aug,...
                        'lam_x0', lam_x_aug,...
                        'p', nlp.p.mpcc_p().val);
                    res_out = tnlp_results;
                    if tnlp_solver.stats.success
                        converged = true;
                    else
                        n_biactive = n_biactive - 1;
                    end 
                    %converged = true;
                end
                w_star = full(tnlp_results.x);
                w_star_orig = w_star(ind_mpcc);

                nu = -full(tnlp_results.lam_x(ind_G));
                xi = -full(tnlp_results.lam_x(ind_H));

                % Index sets
                I_G = find(ind_0p+ind_00);
                I_H = find(ind_p0+ind_00);

                G_new = full(obj.G_fun(w_star_orig, nlp.p.mpcc_p().val));
                H_new = full(obj.H_fun(w_star_orig, nlp.p.mpcc_p().val));

                g_tnlp = full(tnlp_results.g);
                % multipliers (nu, xi already calculated above)
                lam_x_star = full(tnlp_results.lam_x);
                lam_x_star(ind_G) = 0; % Capturing multipliers for non-complementarity 
                lam_x_star(ind_H) = 0; % box constraints
                lam_g_star = tnlp_results.lam_g;
                type_tol = a_tol^2;
            else
                w = nlp.w.mpcc_w(); lbw_orig = nlp.w.mpcc_w().lb; ubw_orig = nlp.w.mpcc_w().ub;
                g = nlp.g.mpcc_g(); lbg = nlp.g.mpcc_g().lb; ubg = nlp.g.mpcc_g().ub;
                p = nlp.p.mpcc_p();
                f = mpcc.f;

                lam_x_aug = lam_x;
                
                g = vertcat(g,10000*G,10000*H); % in case of unlifted complementarity functions scaling the functions by a large value, greatly improves the convergence of the TNLP
                                                % This does not affect the signs of the MPCC multipliers and therefore does not affect the stationarity type verification. 
                lbg = vertcat(lbg,zeros(size(G)),zeros(size(H)));
                ubg = vertcat(ubg,zeros(size(G)),zeros(size(H)));
                lam_g_aug = vertcat(lam_g, zeros(size(G)), zeros(size(H)));

                jac_g = Function('jac_g', {w, p}, {g.jacobian(w)});

                ind_mpcc = 1:length(lam_g);
                ind_G = (1:length(G))+length(lam_g);
                ind_H = (1:length(H))+length(lam_g)+length(G);

                converged = false;
                n_max_biactive = n_biactive;
                n_biactive = n_biactive;
                while ~converged && n_biactive >= 0
                    idx_00 = min_idx(1:n_biactive);
                    ind_00 = false(length(G),1);
                    ind_00(idx_00) = true;
                    
                    ind_0p = G_old<H_old & ~ind_00;
                    ind_p0 = H_old<G_old & ~ind_00;

                    lbG = 0*ones(size(G));
                    ubG = inf*ones(size(G));
                    lbH = 0*ones(size(G));
                    ubH = inf*ones(size(G));
                    ubG(find(ind_00 | ind_0p)) = 1e-12; 
                    ubH(find(ind_00 | ind_p0)) = 1e-12;
                    lbg(ind_G) = lbG;
                    lbg(ind_H) = lbH;
                    ubg(ind_G) = ubG;
                    ubg(ind_H) = ubH;

                    nabla_g = jac_g(w_orig, nlp.p.mpcc_p().val);

                    lbw = lbw_orig;
                    ubw = ubw_orig;
                    w_init = w_orig;

                    
                    casadi_nlp = struct('f', f , 'x', w, 'g', g, 'p', p);
                    opts_casadi_nlp.ipopt.max_iter = 5000;
                    opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
                    opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_bound_frac = 1e-5;
                    opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-5;
                    opts_casadi_nlp.ipopt.tol = obj.opts.complementarity_tol;
                    opts_casadi_nlp.ipopt.acceptable_tol = sqrt(obj.opts.complementarity_tol);
                    opts_casadi_nlp.ipopt.acceptable_dual_inf_tol = sqrt(obj.opts.complementarity_tol);
                    opts_casadi_nlp.ipopt.fixed_variable_treatment = 'relax_bounds';
                    opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
                    opts_casadi_nlp.ipopt.bound_relax_factor = 1e-12;
                    opts_casadi_nlp.ipopt.linear_solver = obj.opts.opts_casadi_nlp.ipopt.linear_solver;
                    opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
                    opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
                    default_tol = 1e-4;
                    opts_casadi_nlp.ipopt.tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
                    opts_casadi_nlp.ipopt.print_level = obj.opts.opts_casadi_nlp.ipopt.print_level;
                    opts_casadi_nlp.print_time = obj.opts.opts_casadi_nlp.print_time;
                    opts_casadi_nlp.ipopt.sb = 'yes';
                    
                    tnlp_solver = nlpsol('tnlp', 'ipopt', casadi_nlp, opts_casadi_nlp);
                    tnlp_results = tnlp_solver('x0', w_init,...
                        'lbx', lbw,...
                        'ubx', ubw,...
                        'lbg', lbg,...
                        'ubg', ubg,...
                        'lam_g0', lam_g_aug,...
                        'lam_x0', lam_x_aug,...
                        'p',nlp.p.mpcc_p().val);
                    res_out = tnlp_results;
                    if tnlp_solver.stats.success
                        converged = true;
                    else
                        n_biactive = n_biactive - 1;
                    end
                    
                    switch tnlp_solver.stats.return_status
                      case {'Solve_Succeeded', 'Solved_To_Acceptable_Level', 'Search_Direction_Becomes_Too_Small'}
                        converged = true;
                      otherwise
                        converged = false;
                    end
                end
                w_star = full(tnlp_results.x);

                nu = -full(tnlp_results.lam_g(ind_G)); % here a - is put, as in casadi Lagrange multipliers for inequality constraints are nonpositive, but we use our defintions with nonnegative multipliers
                xi = -full(tnlp_results.lam_g(ind_H));
                
                % Index sets
                I_G = find(ind_0p+ind_00);
                I_H = find(ind_p0+ind_00);

                G_new = full(obj.G_fun(w_star_orig, nlp.p.mpcc_p().val));
                H_new = full(obj.H_fun(w_star_orig, nlp.p.mpcc_p().val));

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
                    bound = 1.1*(max(abs([nu_biactive;xi_biactive]))+1e-10); % 1e-10 here is added for better plots when all biactive point multipliers are small.  

                    figure()
                    scatter(nu_biactive, xi_biactive, 50, 'o', 'LineWidth', 2);
                    xline(0,'k-.');
                    yline(0,'k-.');
                    xlim([-bound,bound]);
                    ylim([-bound,bound]);
                    xlabel('$\nu$','interpreter','latex')
                    ylabel('$\xi$','interpreter','latex')
                    title('Lagrange multipliers for the comp. constraints')
                    
                    grid on;

                    if all(nu_biactive > -type_tol & xi_biactive > -type_tol)
                        stat_type = "S";
                        if obj.opts.print_level >= 1
                            disp("Converged to S-stationary point.")
                        end
                    elseif all((nu_biactive > -type_tol & xi_biactive > -type_tol) | (abs(nu_biactive.*xi_biactive) < type_tol))
                        stat_type = "M";
                        if obj.opts.print_level >= 1
                            disp("Converged to M-stationary point.")
                        end
                    elseif all(nu_biactive.*xi_biactive > -type_tol)
                        stat_type = "C";
                        if obj.opts.print_level >= 1
                            disp("Converged to C-stationary point.")
                        end
                    elseif all(nu_biactive > -type_tol | xi_biactive > -type_tol)
                        stat_type = "A";
                        if obj.opts.print_level >= 1
                            disp("Converged to A-stationary point.")
                        end
                    else
                        stat_type = "W";
                        if obj.opts.print_level >= 1
                            disp("Converged to W-stationary point, or something has gone wrong.")
                        end
                    end
                else
                    stat_type = "S";
                    if obj.opts.print_level >= 1
                        disp("Converged to S-stationary point.")
                    end
                end
              otherwise
                stat_type = "?";
                if obj.opts.print_level >= 1
                    disp("Could not converge to point from the end of homotopy.");
                end
            end
            
            % output tnlp results
            res_out = tnlp_results;
            w_polished = zeros(size(nlp.w));
            w_polished = w_star_orig;
            res_out.x = w_polished;
        end

        function [solution, improved_point, b_stat] = check_b_stationarity(obj, x0, exitfast)
            import casadi.*
            if ~exist('exitfast', 'var')
                exitfast = true;
            end
            
            mpcc = obj.mpcc;
            nlp = obj.nlp;
            f = mpcc.f;
            x = nlp.w.mpcc_w(); lbx = nlp.w.mpcc_w().lb; ubx = nlp.w.mpcc_w().ub;
            g = nlp.g.mpcc_g(); lbg = nlp.g.mpcc_g().lb; ubg = nlp.g.mpcc_g().ub;
            G = mpcc.G;
            H = mpcc.H;
            G_old = full(obj.G_fun(x0, nlp.p.mpcc_p().val));
            H_old = full(obj.H_fun(x0, nlp.p.mpcc_p().val));

            % We take the square root of the complementarity tolerance for the activity tolerance as it is the upper bound for G, H, for any of the smoothing approaches.
            % i.e. in the worst case G=sqrt(complementarity_tol) and H=sqrt(comp_tol).
            % We add 1e-14 (a value near machine precision) in order to accurately identify points on the boundary.
            if obj.opts.homotopy_steering_strategy == HomotopySteeringStrategy.DIRECT
                a_tol = sqrt(obj.opts.complementarity_tol) + 1e-14;
            else
                a_tol = sqrt(max(nlp.w.s_elastic().res)) + 1e-14;
            end
            ind_00 = G_old<a_tol & H_old<a_tol;
            n_biactive = sum(ind_00);
            if n_biactive == 0 && exitfast
                b_stat = true;
                solution = [];
                improved_point = x0;
                disp("Converged to B-stationary point.")
                return
            end
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
            G_init(ind_00 | ind_0p) = 0;
            H_init = H_old;
            H_init(ind_00 | ind_p0) = 0;
            x0 = vertcat(x0, G_init, H_init);

            g = vertcat(g,g_lift);
            lbg = vertcat(lbg,zeros(size(g_lift)));
            ubg = vertcat(ubg,zeros(size(g_lift)));
            
            p = mpcc.p;
            p0 = nlp.p.mpcc_p().val;

            solver_settings.max_iter = 20;
            solver_settings.tol = 1e-6;
            solver_settings.Delta_TR_init = 10;
            solver_settings.Delta_TR_min = 1e-4;
            solver_settings.tighten_bounds_in_lpcc = false;
            solver_settings.BigM = 1e2;
            if obj.opts.print_level >= 4
                solver_settings.verbose_solver = 1;
            else
                solver_settings.verbose_solver = 0;
            end
            
            if ~n_biactive
                solver_settings.fixed_y_lpcc = y_lpcc;
            end
            %% The problem
            lpec_data = struct('x', x, 'f', f, 'g', g, 'comp1', lift_G, 'comp2', lift_H, 'p', p);
            solver_initalization = struct('x0', x0, 'lbx', lbx, 'ubx', ubx,'lbg', lbg, 'ubg', ubg, 'p', p0, 'y_lpcc', y_lpcc);
            solution = b_stationarity_oracle(lpec_data,solver_initalization,solver_settings);
            
            improved_point = solution.x(ind_mpcc);

            b_stat = ~solution.oracle_status;
            if obj.opts.print_level >= 1
                if b_stat
                    disp("Converged to B-stationary point.")
                else
                    disp("Converged to a point that is not B-stationary.")
                end
            end
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            % TODO(@anton) this returns mpccresults struct from homotopy
            import casadi.*;
            p = inputParser;
            addParameter(p, 'x0', []);
            addParameter(p, 'lbx', []);
            addParameter(p, 'ubx', []);
            addParameter(p, 'lbg', []);
            addParameter(p, 'ubg', []);
            addParameter(p, 'p', []);
            addParameter(p, 'lam_g0', []);
            addParameter(p, 'lam_x0', []);
            % TODO(@anton) Also perhaps initial complementarity multipliers but need to think about
            %              how exactly to implement that.
            parse(p, index_op(1).Indices{:});

            % data and local variables
            nlp = obj.nlp;
            mpcc = obj.mpcc;
            opts = obj.opts;
            plugin = obj.plugin;
            G_fun = Function('G', {mpcc.x, mpcc.p}, {mpcc.G});
            H_fun = Function('H', {mpcc.x, mpcc.p}, {mpcc.H});
            obj.G_fun = G_fun;
            obj.H_fun = H_fun;
            comp_res_fun = Function('comp_res', {mpcc.x, mpcc.p}, {mmax(mpcc.G.*mpcc.H)});
            obj.comp_res_fun = comp_res_fun;
            f_mpcc_fun = Function('f_mpcc', {mpcc.x, mpcc.p}, {mpcc.f});
            obj.f_mpcc_fun = f_mpcc_fun;
            g_mpcc_fun = Function('g_mpcc', {mpcc.x, mpcc.p}, {mpcc.g});
            % Update nlp data
            if ~isempty(p.Results.x0)
                nlp.w.mpcc_w().init = p.Results.x0;
            end
            if ~isempty(p.Results.lbx)
                nlp.w.mpcc_w().lb = p.Results.lbx;
            end
            if ~isempty(p.Results.ubx)
                nlp.w.mpcc_w().ub = p.Results.ubx;
            end
            if ~isempty(p.Results.lbg)
                nlp.g.mpcc_g().lb = p.Results.lbg;
            end
            if ~isempty(p.Results.ubg)
                nlp.g.mpcc_g().ub = p.Results.ubg;
            end
            if ~isempty(p.Results.p)
                nlp.p.mpcc_p().val = p.Results.p;
            end
            if ~isempty(p.Results.lam_g0)
                nlp.g.mpcc_g().init_mult = p.Results.lam_g0;
            end
            if ~isempty(p.Results.lam_x0)
                nlp.w.mpcc_w().init_mult = p.Results.lam_x0;
            end

            if ~opts.assume_lower_bounds % Lower bounds on G, H, not already present in MPCC
                lb = nlp.w.mpcc_w().lb;
                lb(obj.ind_map_G) = 0;
                lb(obj.ind_map_H) = 0;
                nlp.w.mpcc_w().lb = lb;
            end
 

            % Initial conditions
            sigma_k = opts.sigma_0;

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.wall_time = [];
            stats.cpu_time_total = 0;
            stats.wall_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            % TODO(@anton) calculate objective and complementarity residual
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res_fun(nlp.w.mpcc_w().init, nlp.p.mpcc_p().val))];

            % Initialize Results struct
            mpcc_results = struct;
            mpcc_results.W = nlp.w.mpcc_w().init;
            mpcc_results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;
            last_iter_failed = 0;
            timeout = 0;

            if opts.print_level >= 3
                plugin.print_nlp_iter_header();
            end

            while (abs(complementarity_iter) > opts.complementarity_tol || last_iter_failed) &&...
                    ii < opts.N_homotopy &&...
                    (sigma_k > opts.sigma_N || ii == 0) &&...
                    ~timeout
                % homotopy parameter update
                last_iter_failed = 0;
                if ii == 0
                    sigma_k = opts.sigma_0;
                else
                    if isequal(opts.homotopy_update_rule,'linear')
                        sigma_k = opts.homotopy_update_slope*sigma_k;
                    elseif isequal(opts.homotopy_update_rule,'superlinear')
                        sigma_k = max(opts.sigma_N,min(opts.homotopy_update_slope*sigma_k,sigma_k^opts.homotopy_update_exponent));
                    else
                        % TODO(@anton) make this not necessary
                        error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
                    end
                end
                stats.sigma_k = [stats.sigma_k, sigma_k];
                nlp.p.sigma_p().val = sigma_k;

                if ii ~= 0
                    % TODO(@anton) Lets push on casadi devs to allow for changing of options after construction
                    %              (or maybe do it ourselves) and then remove this hack
                    if obj.opts.timeout_cpu
                        plugin.construct_solver(nlp, opts, opts.timeout_cpu - stats.cpu_time_total);
                    elseif obj.opts.timeout_wall
                        plugin.construct_solver(nlp, opts, opts.timeout_wall - stats.wall_time_total);
                    end % HACK ENDS HERE
                end
                start = cputime;
                tic;
                [nlpsol_stats, nlpsol_results] = nlp.solve();
                cpu_time_iter = cputime - start;
                wall_time_iter = toc;
                nlp.w.init = nlp.w.res;
                
                solver_stats = plugin.cleanup_solver_stats(nlpsol_stats);
                
                stats.solver_stats = [stats.solver_stats, solver_stats];

                mpcc_results.nlp_results = [mpcc_results.nlp_results, nlpsol_results];

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

                if opts.timeout_cpu && (stats.cpu_time_total > opts.timeout_cpu)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                if opts.timeout_wall && (stats.wall_time_total > opts.timeout_wall)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                % update results output.
                mpcc_results.W = [mpcc_results.W,nlp.w.mpcc_w().res]; % all homotopy iterations

                % update complementarity and objective stats
                complementarity_iter = full(comp_res_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(f_mpcc_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
                stats.objective = [stats.objective, objective];

                % update counter
                ii = ii+1;
                % Verbose
                if opts.print_level >= 3
                    plugin.print_nlp_iter_info(stats)
                end
            end
            
            mpcc_results.f = full(f_mpcc_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
            mpcc_results.x = nlp.w.mpcc_w().res;
            mpcc_results.p = nlp.p.mpcc_p().val;
            mpcc_results.lam_x = nlp.w.mpcc_w().mult;
            mpcc_results.g = nlp.g.mpcc_g().eval;
            mpcc_results.lam_g = nlp.g.mpcc_g().mult;
            % TODO(@anton) it may not always be possible to calculate lam_G and lam_H
            %              figure out when this is possible.
            mpcc_results.G = full(G_fun(mpcc_results.x, mpcc_results.p));
            % mpcc_results.lam_G = ;
            mpcc_results.H = full(H_fun(mpcc_results.x, mpcc_results.p));
            % mpcc_results.lam_H = ;
            
            stats.converged = obj.complementarity_tol_met(stats) && ~last_iter_failed && ~timeout;
            stats.success = stats.converged;
            stats.constraint_violation = obj.compute_constraint_violation(mpcc_results.x, mpcc_results.g);
            
            if opts.calculate_stationarity_type || opts.polishing_step
                if last_iter_failed || ~obj.complementarity_tol_met(stats) || timeout
                    stat_type = "?";
                    disp("Not checking stationarity due to failure of homotopy to converge.");
                else                    
                    [w_polished, res_out, stat_type, n_biactive] = obj.check_multiplier_based_stationarity(nlp.w.mpcc_w().res);
                    [sol, w_polished, b_stat] = obj.check_b_stationarity(w_polished);
                    if stat_type ~= "?"
                        mpcc_results.x = w_polished;
                        mpcc_results.f = full(f_mpcc_fun(w_polished, nlp.p.mpcc_p().val));
                        mpcc_results.g = full(g_mpcc_fun(w_polished, nlp.p.mpcc_p().val));
                        mpcc_results.G = full(G_fun(w_polished, nlp.p.mpcc_p().val));
                        mpcc_results.H = full(H_fun(w_polished, nlp.p.mpcc_p().val));
                        % TODO (@anton) also recalculate multipliers and g?
                        mpcc_results.nlp_results = [mpcc_results.nlp_results, res_out];
                        complementarity_iter = full(obj.comp_res_fun(w_polished, nlp.p.mpcc_p().val));
                        stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                    end
                    stats.stat_type = stat_type;
                    stats.b_stationary = b_stat;
                end
            end            
            varargout{1} = mpcc_results;
            obj.stats = stats;
        end
        
        function obj = parenAssign(obj,index_op,varargin)
            error('Invalid operation');
        end
        
        function obj = parenDelete(obj,index_op)
            error('Invalid operation')
        end

        function n = parenListLength(obj,index_op,ctx)
            n = 1;
        end

        function print_infeasibility(obj, results)
            if obj.opts.print_details_if_infeasible
                print_problem_details(results,obj.model,obj.mpcc, []);
            end
            if obj.opts.pause_homotopy_solver_if_infeasible
                keyboard
            end
        end
    end
end
