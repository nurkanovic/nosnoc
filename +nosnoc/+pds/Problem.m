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
classdef Problem < vdx.problems.Mpcc
    properties (Access=public)
        model
        opts
    end

    methods (Access=public)
        function obj = Problem(model, opts)
            obj = obj@vdx.problems.Mpcc();
            obj.model = model;
            obj.opts = opts;
        end

        function create_variables(obj)
            model = obj.model;
            opts = obj.opts;
            opts.preprocess();
            model.verify_and_backfill(opts);
            model.generate_variables();
            model.generate_functions();
            dims = model.dims;
            
            obj.p.sigma(1) = {{'sigma',1},0,inf,0};
            obj.p.rho_h_p(1) = {{'rho_h_p',1},0,inf,1};
            obj.p.gamma_h(1) = {{'gamma_h',1},0,inf,1e-1};
            obj.p.T(1) = {{'T',1},0,inf,opts.T};
            obj.p.p_global(1) = {model.p_global,-inf, inf, model.p_global_val};

            % other derived values
            t_stage = opts.T/opts.N_stages;
            h0 = opts.h;
            
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}};
            obj.w.lambda(0,0,opts.n_s) = {{['lambda_0'], dims.n_c},0,0};
            obj.w.v_global(1) = {{'v_global',dims.n_v_global}, model.lbv_global, model.ubv_global, model.v0_global};
            for ii=1:opts.N_stages
                obj.w.u(ii) = {{['u_' num2str(ii)], dims.n_u}, model.lbu, model.ubu, model.u0};
                obj.p.p_time_var(ii) = {{['p_time_var_' num2str(ii)], dims.n_p_time_var}, -inf, inf, model.p_time_var_val};
                if obj.opts.use_speed_of_time_variables
                    obj.w.sot(ii) = {{['sot_' num2str(ii)], 1}, 0, 100, 1};
                end
                for jj=1:opts.N_finite_elements
                    if obj.opts.use_fesd
                        obj.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, (1-opts.gamma_h)*h0, (1+opts.gamma_h)*h0, h0};
                    end
                    if (strcmp(obj.opts.step_equilibration,'linear')||...
                        strcmp(obj.opts.step_equilibration,'linear_tanh')||...
                        strcmp(obj.opts.step_equilibration,'linear_relaxed')) && jj > 1
                        obj.w.B_max(ii,jj) ={{['B_max_' num2str(ii) '_' num2str(jj)], dims.n_c},-inf,inf};
                        obj.w.pi_lambda(ii,jj) ={{['pi_lambda_' num2str(ii) '_' num2str(jj)], dims.n_c},-inf,inf};
                        obj.w.pi_c(ii,jj) ={{['pi_c_' num2str(ii) '_' num2str(jj)], dims.n_c},-inf,inf};
                        obj.w.lambda_lambda(ii,jj) ={{['lambda_lambda_' num2str(ii) '_' num2str(jj)], dims.n_c},0,inf};
                        obj.w.lambda_c(ii,jj) ={{['lambda_c_' num2str(ii) '_' num2str(jj)], dims.n_c},0,inf};
                        obj.w.eta(ii,jj) ={{['eta_' num2str(ii) '_' num2str(jj)], dims.n_c},0,inf};
                        obj.w.nu(ii,jj) ={{['nu_' num2str(ii) '_' num2str(jj)], 1},0,inf};
                    end
                    for kk=1:opts.n_s
                        obj.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], dims.n_x}, model.lbx, model.ubx, model.x0};
                        obj.w.z(ii,jj,kk) = {{['z_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], dims.n_z}, model.lbz, model.ubz, model.z0};
                        obj.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], dims.n_c},0,inf};
                    end
                end
            end
        end

        function forward_sim_constraints(obj)
            import casadi.*
            model = obj.model;
            opts = obj.opts;
            if obj.opts.use_fesd
                t_stage = obj.p.T(1)/opts.N_stages;
                h0 = obj.p.T(1).init/(opts.N_stages*opts.N_finite_elements);
            else
                h0 = obj.p.T(1).init/(opts.N_stages*opts.N_finite_elements);
            end
            
            % Define functions from obj.data
            nabla_c = model.c.jacobian(model.x)';
            E = model.E;
            v_global = obj.w.v_global(1);
            p_global = obj.p.p_global(1);
            c_fun = model.c_fun;

            x_prev = obj.w.x(0,0,opts.n_s);
            for ii=1:opts.N_stages
                ui = obj.w.u(ii);
                p_stage = obj.p.p_time_var(ii);
                p =[p_global;p_stage];
                if obj.opts.use_speed_of_time_variables
                    s_sot = obj.w.sot(ii);
                else
                    s_sot = 1;
                end
                sum_h = 0;
                for jj=1:opts.N_finite_elements
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    else
                        h = h0;
                    end
                    for kk=1:opts.n_s
                        x_ijk = obj.w.x(ii,jj,kk);
                        z_ijk = obj.w.z(ii,jj,kk);
                        lambda_ijk = obj.w.lambda(ii,jj,kk);
                        fj = s_sot*model.f_x_fun(x_ijk, z_ijk, ui, lambda_ijk, v_global, p);
                        qj = s_sot*model.f_q_fun(x_ijk, z_ijk, ui, v_global, p);
                        xk = opts.C_irk(1, kk+1) * x_prev;
                        for rr=1:opts.n_s
                            x_ijr = obj.w.x(ii,jj,rr);
                            xk = xk + opts.C_irk(rr+1, kk+1) * x_ijr;
                        end
                        obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                        obj.g.algebraics(ii,jj,kk) = {model.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                        % also add non-negativity constraint on c
                        obj.g.c_nonnegative(ii,jj,kk) = {c_fun(x_ijk, v_global, p), 0, inf};
                        %add path constraints
                        obj.g.g_path(ii,jj,kk) = {model.g_path_fun(x_ijk, z_ijk, ui, v_global, p), model.lbg_path, model.ubg_path};
                        % also integrate the objective
                        obj.f = obj.f + opts.B_irk(kk+1)*h*qj;
                    end
                    x_prev = obj.w.x(ii,jj,opts.n_s);
                end
                if obj.opts.use_fesd
                    obj.g.sum_h(ii) = {t_stage-sum_h};
                end
                if obj.opts.use_speed_of_time_variables
                    x_end = obj.w.x(ii,opts.N_finite_elements,opts.n_s);
                    x_start = obj.w.x(0,0,opts.n_s);
                    obj.g.g_equidistant_grid(ii) = {(x_end(end)-x_start(end)) - t_stage*ii};
                end
            end

            % Terminal cost
            obj.f = obj.f + model.f_q_T_fun(obj.w.x(ii,jj,kk), obj.w.z(ii,jj,kk), v_global, p_global);

            % Terminal constraint
            obj.g.terminal(0) = {model.g_terminal_fun(obj.w.x(ii,jj,kk), obj.w.z(ii,jj,kk), v_global, p_global), model.lbg_terminal, model.ubg_terminal};
        end

        function generate_complementarities(obj)
            import casadi.*
            opts = obj.opts;
            model = obj.model;
            % Do Cross-Complementarity
            
            x_prev = obj.w.x(0,0,opts.n_s);
            lambda_prev = obj.w.lambda(0,0,opts.n_s);
            v_global = obj.w.v_global(1);
            p_global = obj.p.p_global(1);
            c_fun = model.c_fun;

            G = [];
            H = [];
            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p =[p_global;p_stage];
                        for jj=1:opts.N_finite_elements
                            Gij = [];
                            Hij = [];
                            for kk=1:opts.n_s
                                for rr=1:opts.n_s
                                    x_ijk = obj.w.x(ii,jj,kk);
                                    lambda_ijr = obj.w.lambda(ii,jj,rr);
                                    Gij = vertcat(Gij,c_fun(x_ijk, v_global, p));
                                    Hij = vertcat(Hij,lambda_ijr);
                                end
                            end
                            for kk=1:opts.n_s 
                                x_ijk = obj.w.x(ii,jj,kk);
                                lambda_ijr = lambda_prev;
                                Gij = vertcat(Gij,c_fun(x_ijk, v_global, p));
                                Hij = vertcat(Hij,lambda_prev);
                            end
                            for rr=1:opts.n_s
                                lambda_ijr = obj.w.lambda(ii,jj,rr);
                                Gij = vertcat(Gij,c_fun(x_prev, v_global, p));
                                Hij = vertcat(Hij,lambda_ijr);
                            end
                            
                            obj.G.cross_comp(ii,jj) = {Gij};
                            obj.H.cross_comp(ii,jj) = {Hij};
                            x_prev = obj.w.x(ii,jj,opts.n_s);
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.FE_STAGE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p =[p_global;p_stage];
                        for jj=1:opts.N_finite_elements
                            Gij = [];
                            Hij = [];
                            sum_c = c_fun(x_prev);
                            for kk=1:opts.n_s
                                x_ijk = obj.w.x(ii,jj,kk);
                                sum_c = sum_c + c_fun(x_ijk, v_global, p);
                            end
                            Gij = vertcat(Gij,sum_c);
                            Hij = vertcat(Hij,lambda_prev);
                            for rr=1:opts.n_s
                                lambda_ijr = obj.w.lambda(ii,jj,rr);
                                Gij = vertcat(Gij,sum_c);
                                Hij = vertcat(Hij,lambda_ijr);
                            end
                            
                            obj.G.cross_comp(ii,jj) = {Gij};
                            obj.H.cross_comp(ii,jj) = {Hij};
                            x_prev = obj.w.x(ii,jj,opts.n_s);
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p =[p_global;p_stage];
                        for jj=1:opts.N_finite_elements
                            Gij = [];
                            Hij = [];
                            sum_lambda = lambda_prev;
                            for rr=1:opts.n_s
                                lambda_ijr = obj.w.lambda(ii,jj,rr);
                                sum_lambda = sum_lambda + lambda_ijr;
                            end
                            Gij = vertcat(Gij,c_fun(x_prev, v_global, p));
                            Hij = vertcat(Hij,sum_lambda);
                            for kk=1:opts.n_s
                                x_ijk = obj.w.x(ii,jj,kk);
                                Gij = vertcat(Gij,c_fun(x_ijk, v_global, p));
                                Hij = vertcat(Hij,sum_lambda);
                            end
                            
                            obj.G.cross_comp(ii,jj) = {Gij};
                            obj.H.cross_comp(ii,jj) = {Hij};
                            x_prev = obj.w.x(ii,jj,opts.n_s);
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.FE_FE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p =[p_global;p_stage];
                        for jj=1:opts.N_finite_elements
                            Gij = c_fun(x_prev, v_global, p);
                            Hij = lambda_prev;
                            for kk=1:opts.n_s
                                x_ijk = obj.w.x(ii,jj,kk);
                                lambda_ijk = obj.w.lambda(ii,jj,kk);
                                Gij = Gij + c_fun(x_ijk, v_global, p);
                                Hij = Hij + lambda_ijk;
                            end
                            obj.G.cross_comp(ii,jj) = {Gij};
                            obj.H.cross_comp(ii,jj) = {Hij};
                            x_prev = obj.w.x(ii,jj,opts.n_s);
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                end
            else
                Gij = [];
                Hij = [];
                for kk=1:opts.n_s
                    x_ijk = obj.w.x(ii,jj,kk);
                    lambda_ijk = obj.w.lambda(ii,jj,kk);
                    Gij = [Gij;c_fun(x_ijk)];
                    Hij = [Hij;lambda_ijk];
                end
                obj.G.comp(ii,jj) = {Gij};
                obj.H.comp(ii,jj) = {Hij};
            end
        end

        function step_equilibration(obj)
            model = obj.model;
            opts = obj.opts;
            h0 = opts.h;
            v_global = obj.w.v_global(1);
            p_global = obj.p.p_global(1);
            c_fun = model.c_fun;
            switch obj.opts.step_equilibration
              case 'heuristic_mean'
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements
                        obj.f = obj.f + obj.p.gamma_h(1)*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case 'heuristic_diff'
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements
                        obj.f = obj.f + obj.p.gamma_h(1)*(obj.w.h(ii,jj)-obj.w.h(ii,jj-1))^2;
                    end
                end
              case 'l2_relaxed_scaled'
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:opts.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk), v_global, p);
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:opts.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk), v_global, p);
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p(1) * tanh(eta/opts.step_equilibration_sigma) * delta_h.^2;
                    end
                end
              case 'l2_relaxed'
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:opts.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk), v_global, p);
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:opts.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk), v_global, p);
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p(1) * eta * delta_h.^2
                    end
                end
              case 'direct'
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:opts.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk), v_global, p);
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:opts.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk), v_global, p);
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.g.step_equilibration(ii,jj) = {eta*delta_h, 0, 0};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
              case 'direct_homotopy'
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:opts.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk), v_global, p);
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:opts.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk), v_global, p);
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        homotopy_eq = [eta*delta_h - sigma;eta*delta_h + sigma];
                        obj.g.step_equilibration(ii,jj) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
              case 'mlcp'
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements
                        sigma_c_b = 0;
                        sigma_lambda_b = 0;
                        for kk=1:opts.n_s
                            sigma_c_b = sigma_c_b + c_fun(obj.w.x(ii,jj-1,kk), v_global, p);
                            sigma_lambda_b = sigma_lambda_b + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_f = 0;
                        sigma_lambda_f = 0;
                        for kk=1:opts.n_s
                            sigma_c_f = sigma_c_f + c_fun(obj.w.x(ii,jj,kk), v_global, p);
                            sigma_lambda_f = sigma_lambda_f + obj.w.lambda(ii,jj,kk);
                        end

                        
                        % todo ideally we output G and H instead of doing all of the stuff here but ok.
                        lambda_c = obj.w.lambda_c(ii,jj);
                        lambda_lambda = obj.w.lambda_lambda(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_f;pi_c-sigma_c_b;sigma_c_f+sigma_c_b-pi_c],0,inf};
                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_f;pi_lambda-sigma_lambda_b;sigma_lambda_f+sigma_lambda_b-pi_lambda],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_lambda-lambda_c;
                            B_max-pi_c;
                            B_max-pi_lambda;
                            (B_max-pi_c).*lambda_c - sigma;
                            (B_max-pi_lambda).*lambda_lambda - sigma];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1);-inf*ones(n_c,1);-inf*ones(n_c,1)],
                            [0*ones(n_c,1);inf*ones(n_c,1);inf*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1)]};

                        % eta calculation
                        eta_const = [eta-pi_lambda;eta-pi_c;eta-pi_lambda-pi_c+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(n_c,1);-inf*ones(n_c,1);zeros(n_c,1)],
                            [zeros(n_c,1);zeros(n_c,1);inf*ones(n_c,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        %M = 1e5;
                        M=t_stage;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + (1/h0)*nu*M;
                            delta_h - (1/h0)*nu*M];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                    end
                end
            end
        end

        function create_solver(obj, casadi_options)
            obj.forward_sim_constraint();
            obj.generate_complementarities();
            obj.step_equilibration();
            create_solver@vdx.problems.Mpcc(obj);
        end
    end
end
