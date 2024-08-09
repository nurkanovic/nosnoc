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
clear all;

import casadi.*

F_friction = 2;

%% model stuff
model = get_cart_pole_with_friction_model(0, F_friction);

%%
n_s = 1;
N = 30; % number of control intervals
T = 5; % time horizon
N_FE = 2; % number of integration steps = Finite elements, without switch detections
x_ref = [0; 180/180*pi; 0; 0]; % target position

%% NLP formulation & solver creation
nlp = setup_collocation_nlp(model, T, N, n_s, N_FE);
casadi_nlp = struct('f', nlp.f, 'x', nlp.w, 'p', nlp.p, 'g', nlp.g);
solver = nlpsol('solver', 'ipopt', casadi_nlp);

%% solve nlp
sigma0 = 1e-0;
% Solve the NLP
sol = solver('x0', nlp.w0, 'lbx', nlp.lbw, 'ubx', nlp.ubw,...
    'lbg', nlp.lbg, 'ubg', nlp.ubg, 'p', sigma0);
w_opt = full(sol.x);

%% extract results
n_x = length(model.x);
n_u = length(model.u);

idx_diff = n_u + n_x * (n_s * N_FE) + (N_FE)*n_x; % between x values at shooting nodes
q1_opt = w_opt(1:idx_diff:end)';
q2_opt = w_opt(2:idx_diff:end)';
v1_opt = w_opt(3:idx_diff:end)';
v2_opt = w_opt(4:idx_diff:end)';
u_opt = w_opt(5:idx_diff:end)';

results = struct();
results.x = [q1_opt; q2_opt; v1_opt; v2_opt];
results.t_grid = linspace(0, T, N+1);
results.t_grid_u = linspace(0, T, N+1);
results.u = u_opt;

plot_cart_pole_smoothing_trajectory(results, T/N, x_ref);

%%
figure;
hold on
F_friction_opt = F_friction*tanh(v1_opt/sigma0);
plot(results.t_grid, F_friction_opt)
stairs(results.t_grid, [u_opt, nan]');

grid on
yline(0, 'k--');
ylabel('$F_{\mathrm{fr}}(t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
legend('smoothed friction forces', 'controls')

%% objective experiment
if 0
warm_starting = 0;

n_sigmas = 15;
sigma_values = logspace(1, -8, n_sigmas);
f_opt = [];
w0 = nlp.w0;
for sigma = sigma_values
    sol = solver('x0', w0, 'lbx', nlp.lbw, 'ubx', nlp.ubw,...
        'lbg', nlp.lbg, 'ubg', nlp.ubg, 'p', sigma);
    f_opt = [f_opt, full(sol.f)];
    if warm_starting
        w0 = full(sol.x);
    end
end

% plot
figure
semilogx(sigma_values, f_opt)
grid on
xlabel('$\sigma$ smoothing', 'Interpreter', 'latex')
ylabel('objective', 'Interpreter', 'latex')
end

function nlp = setup_collocation_nlp(model, T, N, n_s, N_FE)
    import casadi.*
    [B, C, D] = generate_collocation_coefficients(n_s);

    %% Extract from model
    L = model.f_q;
    n_x = length(model.x);
    n_u = length(model.u);

    % Continuous time dynamics
    f = Function('f', {model.x, model.u, model.p}, {model.f_expl_ode, L});
    f_terminal = Function('f', {model.x,}, {model.f_q_T});

    % Control discretization
    h = T/(N*N_FE);
    x0 = model.x0;

    %% Casadi NLP formulation
    w = {};
    w0 = [];
    lbw = [];
    ubw = [];
    objective = 0;
    g = {};
    lbg = [];
    ubg = [];

    % "Lift" initial conditions
    Xk = SX.sym('X0', n_x);
    w = {w{:}, Xk};
    lbw = [lbw; x0];
    ubw = [ubw; x0];
    w0 = [w0; x0];

    % Loop over shooting intervals
    for k=0:N-1
        % New NLP variable for the control
        Uk = SX.sym(['U_' num2str(k)],n_u);
        w = {w{:}, Uk};
        lbw = [lbw; model.lbu];
        ubw = [ubw; model.ubu];
        w0 = [w0;  0];

        % Loop over integration steps / finite elements
        for i_fe=1:N_FE
            Xk_end = D(1)*Xk;

            % State at collocation points
            Xki = {};
            for i=1:n_s
                Xki{i} = SX.sym(['X_' num2str(k) '_' num2str(i_fe) '_' num2str(i)], n_x);
                w = {w{:}, Xki{i}};
                lbw = [lbw; model.lbx];
                ubw = [ubw; model.ubx];
                w0 = [w0; x0];
            end
            % Loop over collocation points
            for i=1:n_s
                % Expression for the state derivative at the collocation point
                xp = C(1,i+1)*Xk;
                for r=1:n_s
                    xp = xp + C(r+1,i+1)*Xki{r};
                end

                % Append collocation equations
                [fi, qi] = f(Xki{i}, Uk, model.p);
                g = {g{:}, h*fi - xp};
                lbg = [lbg; zeros(n_x,1)];
                ubg = [ubg; zeros(n_x,1)];

                % Add contribution to the end state
                Xk_end = Xk_end + D(i+1)*Xki{i};

                % Add contribution to quadrature function
                objective = objective + B(i+1)*qi*h;
            end

            % New NLP variable for state at end of interval
            Xk = SX.sym(['X_' num2str(k+1)], n_x);
            w = {w{:}, Xk};
            lbw = [lbw; model.lbx];
            ubw = [ubw; model.ubx];
            w0 = [w0; x0];

            % Add equality constraint
            g = {g{:}, Xk_end-Xk};
            lbg = [lbg; zeros(n_x,1)];
            ubg = [ubg; zeros(n_x,1)];
        end
    end

    objective = objective + f_terminal(Xk);

    nlp = struct();
    nlp.f = objective;
    nlp.w = vertcat(w{:});
    nlp.w0 = w0;
    nlp.lbw = lbw;
    nlp.ubw = ubw;
    nlp.p = model.p;
    nlp.g = vertcat(g{:});
    nlp.lbg = lbg;
    nlp.ubg = ubg;
end


function [B, C, D] = generate_collocation_coefficients(n_s)
    import casadi.*
    %% collocation using casadi
    % Get collocation points
    tau_root = [0 collocation_points(n_s, 'radau')];

    % Coefficients of the collocation equation
    C = zeros(n_s+1,n_s+1);

    % Coefficients of the continuity equation
    D = zeros(n_s+1, 1);

    % Coefficients of the quadrature function
    B = zeros(n_s+1, 1);

    % Construct polynomial basis
    for j=1:n_s+1
        % Construct Lagrange polynomials to get the polynomial basis at the collocation point
        coeff = 1;
        for r=1:n_s+1
            if r ~= j
                coeff = conv(coeff, [1, -tau_root(r)]);
                coeff = coeff / (tau_root(j)-tau_root(r));
            end
        end
        % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        D(j) = polyval(coeff, 1.0);

        % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        pder = polyder(coeff);
        for r=1:n_s+1
            C(j,r) = polyval(pder, tau_root(r));
        end

        % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        pint = polyint(coeff);
        B(j) = polyval(pint, 1.0);
    end
end
