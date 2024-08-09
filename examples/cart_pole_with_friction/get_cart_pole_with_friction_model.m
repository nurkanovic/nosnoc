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

function model = get_cart_pole_with_friction_model(use_nosnoc, F_friction)
    import casadi.*
    if use_nosnoc
        model = nosnoc.model.Pss();
    else
        model = struct();
    end

    % fixed values
    m1 = 1; % cart
    m2 = 0.1; % link
    g = 9.81;
    link_length = 1;

    % symbolics
    px = SX.sym('px');
    theta = SX.sym('theta');
    q = vertcat(px, theta);

    v = SX.sym('v');
    theta_dot = SX.sym('theta_dot');
    q_dot = vertcat(v, theta_dot);

    x = vertcat(q, q_dot);

    u = SX.sym('u'); % control

    % Inertia matrix
    M = [m1 + m2, m2*link_length*cos(theta);...
        m2 *link_length*cos(theta),  m2*link_length^2];
    % Coriolis force
    C = [0, -m2 * link_length*theta_dot*sin(theta);...
        0,   0];

    % all forces = Gravity + Control + Coriolis (+Friction)
    f_all = [0; -m2*g*link_length*sin(theta)] + [u; 0] - C*q_dot;

    if use_nosnoc
        % switching function c: cart velocity
        model.c = v;
        % Sign matrix % f_1 for c>0, f_2 for c<0
        model.S = [1; -1];

        % Dynamics for v > 0
        f_1 = [q_dot;...
                inv(M)*(f_all-[F_friction; 0])];
        % Dynamics for v<0
        f_2 = [q_dot;...
                inv(M)*(f_all+[F_friction; 0])];
        model.F = [f_1, f_2];
    else
        sigma = SX.sym('sigma');
        model.p = sigma;
        model.f_expl_ode = [q_dot; inv(M) * (f_all - [F_friction*tanh(v/ sigma); 0])];
    end

    % specify initial and desired state
    x0 = [1; 0/180*pi; 0; 0]; % start downwards
    x_ref = [0; 180/180*pi; 0; 0]; % end upwards

    % bounds
    ubx = [5; inf; inf; inf];
    lbx = [-5; -inf; -inf; -inf];
    u_max = 30;

    % cost
    Q = diag([0.01;0;0;0]);%diag([10; 100; 1; 1]);
    R = 1;
    model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;

    % terminal cost could be added
    Q_terminal = diag([500; 5000; 200; 200]);
    model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);

    model.lbx = lbx;
    model.ubx = ubx;
    model.x = x;
    model.x0 =  x0;
    model.u = u;
    model.lbu = -u_max;
    model.ubu = u_max;
end
