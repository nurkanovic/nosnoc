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

%
%
function [varargout] = homotopy_solver(varargin)
% homotopy_solver(solver,model,settings,solver_initialization)
solver = varargin{1};
model = varargin{2};
settings = varargin{3};
solver_initialization = varargin{4};

if nargin>4
    model_int = varargin{5};
    settings_int = varargin{6};
end

import casadi.*
%%  unfold data
settings_bkp = settings;
unfold_struct(settings,'caller')
settings = settings_bkp;
unfold_struct(solver_initialization,'caller')

comp_res = model.comp_res;
nabla_J_fun = model.nabla_J_fun;
s_elastic_iter = 1;

sigma_k = sigma_0;
x0 = solver_initialization.w0(1:model.dimensions.n_x);

% cont algebraics initialization
lambda00 = [];
gamma_00 = [];
p_vt_00 = [];
n_vt_00  = [];
gamma_d00 = [];
delta_d00 = [];
y_gap00 = [];
switch settings.dcs_mode
    case 'Stewart'
        lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
    case 'Step'
        lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
    case 'CLS'
        % TODO: reconsider this if 0th element has an impulse
        y_gap00 = f_c_fun(x0);
        if model.friction_exists
            switch settings.friction_model
                case 'Polyhedral'
                    v0 = x0(n_q+1:end);
                    D_tangent_0 = D_tangent_fun(x0);
                    v_t0 = D_tangent_0'*v0;
                    for ii = 1:model.dims.n_contacts
                        ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                        gamma_d00 = [gamma_d00;norm(v_t0(ind_temp))/model.dims.n_t];
                        delta_d00 = [delta_d00;D_tangent_0(:,ind_temp)'*v_t0+gamma_d00(ii)];
                    end
                case 'Conic'
                    v0 = x0(n_q+1:end);
                    v_t0 = J_tangent_fun(x0)'*v0;
                    for ii = 1:n_contacts
                        ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                        v_ti0 = v_0(ind_temp);
                        gamma_00 = [gamma_00;norm(v_ti0)];
                        p_vt_00 = [p_vt_00;max(v_ti0,0)];
                        n_vt_00 = [n_vt_00;max(-v_ti0,0)];
                    end
            end
        end
end
p_val = [model.p_val(:);x0(:);lambda00(:);y_gap00(:);gamma_00(:);gamma_d00(:);p_vt_00(:);n_vt_00(:)];

complementarity_stats = [full(comp_res(w0, p_val))];

cpu_time = [];
homotopy_iterations = [];
w0_base = w0;
W = [w0];

lbw_h = lbw; ubw_h = ubw;
lbw_h(model.ind_h) = model.h_k(1);
ubw_h(model.ind_h) = model.h_k(1);

%% homotopy loop
complementarity_iter = 1;
ii = 0;

if print_level >= 3
    fprintf('\niter\t\tsigma\t\tcompl_res\tobjective\tCPU time\tNLP iters\tstatus\n')
end

while (complementarity_iter) > comp_tol && ii < N_homotopy && (sigma_k > sigma_N || ii == 0)
    % homotopy parameter update
    if ii == 0
        sigma_k = sigma_0;
    else
        if isequal(homotopy_update_rule,'linear')
            sigma_k = homotopy_update_slope*sigma_k;
        elseif isequal(homotopy_update_rule,'superlinear')
            sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent));
        else
            error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
        end
    end
    p_val(1) = sigma_k;
    if h_fixed_to_free_homotopy
        p_val(3) = 1+(sigma_k*1e4);
    end

    % solve problem with fixed step size
    if h_fixed_iterations && use_fesd  && ii < h_fixed_max_iter
        tic
        results = solver('x0', w0, 'lbx', lbw_h, 'ubx', ubw_h,'lbg', lbg, 'ubg', ubg,'p',p_val);
        cpu_time_iter = toc ;
        w_opt = full(results.x);
        if ~h_fixed_change_sigma
            ii = -1; h_fixed_iterations  = 0;
        end
    else
        tic
        results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
        cpu_time_iter = toc ;
    end
    if isequal(solver.stats.return_status,'Infeasible_Problem_Detected')
        warning('nosnoc:homotopy_solver:NLP_infeasible', 'NLP infeasible: try different mpcc_mode or check problem functions.');
    end

    cpu_time = [cpu_time,cpu_time_iter];
    w_opt = full(results.x);
    w0 = w_opt;
    W = [W,w_opt]; % all homotopy iterations

    % complementarity
    complementarity_iter = full(comp_res(w_opt, p_val));
    complementarity_stats = [complementarity_stats;complementarity_iter];

    objective = full(model.problem.objective_fun(w_opt, p_val));
    % update counter
    ii = ii+1;

    % Verbose
    if print_level >= 3
        fprintf('%d\t\t%2.2e\t%2.2e\t%.3f\t\t%.3f\t\t%d\t\t%s\n',ii, sigma_k, complementarity_iter, objective,...
            cpu_time_iter, solver.stats.iter_count, solver.stats.return_status);
    end
    %
    %     if complementarity_iter> 1e1 && ii >= ratio_for_homotopy_stop*N_homotopy
    %         error('The homotopy loop is diverging. Try chaning parameters of the homotopy loop or check is the OCP well posed.')
    %         break;
    %     end
end

%% polish homotopy solution with fixed active set.
if polishing_step
    [results] = polishing_homotopy_solution(model,settings,results,sigma_k);
    %     [results] = polishing_homotopy_solution(model,settings,results,sigma_k,solver,solver_initialization);
    complementarity_iter = results.complementarity_iter;
    complementarity_stats = [complementarity_stats;complementarity_iter];
    W = [W,results.w_opt];
end

%% collect stats
results.W = W;
stats.complementarity_stats = complementarity_stats;
stats.cpu_time = cpu_time;
stats.cpu_time_total = sum(cpu_time);
stats.sigma_k = sigma_k;
stats.homotopy_iterations = ii;

%% loop output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = solver_initialization;
end

