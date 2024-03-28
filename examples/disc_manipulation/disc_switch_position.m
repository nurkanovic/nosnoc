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

%% Switch positions
 
clear all;
close all;
clc;
import casadi.*

filename = 'discs_switch_position.gif';
%%
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 1;  
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 1;
problem_options.cross_comp_mode = 7;
solver_options.mpcc_mode = MpccMode.Scholtes_ineq; 
solver_options.homotopy_update_slope = 0.25; 
solver_options.N_homotopy = 100;
solver_options.comp_tol = 1e-6;
solver_options.opts_casadi_nlp.ipopt.max_iter = 5e3;
solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-5;
% problem_options.g_path_at_fe = 1;
% problem_options.g_path_at_stg = 1;
solver_options.print_level = 5;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
%solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretizatioon
T = 3;
N_stg = 20; % control intervals
N_FE = 3;  % integration steps per control interval

%% model parameters
m1 = 2;
m2 = 1;
r1 = 0.3;
r2 = 0.2;

       
q10 = [-1; 0];
q20 = [1;0];
v10 = [0;0];
v20 = [0;0];

q_target1 = q20;
q_target2 = q10;

x0 = [q10;q20;v10;v20];
ubx = [10; 10;10; 10; 5; 5; 5; 5]; 
lbx = -ubx;
ubu = [20;20];
lbu= -ubu;

x_ref = [q_target1;q_target2;zeros(4,1)];
u_ref = [0;0];

Q = diag([5;5;10;10;0*ones(4,1)]);
R = diag([0.1 0.1]);
Q_terminal = 100*Q;
 
%% Symbolic variables and bounds
q = SX.sym('q',4);
v = SX.sym('v',4); 
u = SX.sym('u',2);

q1 = q(1:2);
q2 = q(3:4);
v1 = v(1:2);
v2 = v(3:4);
x = [q;v];

model = NosnocModel();
problem_options.T = T;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu_f = 0;
model.a_n = 10;
model.x0 = x0; 
model.dims.n_dim_contact = 2;

cv = 2;
eps = 1e-1;
f_drag = cv*[v1/norm(v1+eps);v2/norm(v2+eps)];

model.M = diag([m1;m1;m2;m2]); % inertia/mass matrix;
model.f_v = [u;...
           zeros(2,1)]-f_drag;
model.f_c = [norm(q1-q2)^2-(r1+r2)^2];
% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
%% Objective
model.f_q = 1*(x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
%% Call nosnoc solver
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%% read and plot results
p1 = results.x(1,:);
p2 = results.x(2,:);
p3 = results.x(3,:);
p4 = results.x(4,:);
v1 = results.x(5,:);
v2 = results.x(6,:);
v3 = results.x(7,:);
v4 = results.x(8,:);
t_opt = results.x(9,:);

%% animation
figure('Renderer', 'painters', 'Position', [100 100 1000 800])
x_min = min([p1,p2,p3,p4])-1;
x_max = max([p1,p2,p3,p4])+1;

tt = linspace(0,2*pi,100);
x_t = cos(tt);
y_t = sin(tt);

for ii = 1:length(p1)
    plot(r1*x_t+p1(ii),r1*y_t+p2(ii),'k-','LineWidth',2);
    hold on
    plot(r2*x_t+p3(ii),r2*y_t+p4(ii),'r-','LineWidth',2);
    plot(r1*x_t+q_target1(1),r1*y_t+q_target1(2),'color',[0 0 0 0.6]);
    plot(r2*x_t+q_target2(1),r2*y_t+q_target2(2),'color',[1 0 0 0.6]);
    axis equal
    xlim([x_min x_max])
    ylim([x_min x_max])
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('$y$ [m]','Interpreter','latex');

    % save gif
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',problem_options.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',problem_options.h_k(1));
    end

    if ii~=length(p1)
        clf;
    end
end

%%

if 1
%     figure('Renderer', 'painters', 'Position', [100 100 1400 600])
    figure
    subplot(311)
    plot(t_opt,v1,'LineWidth',1.5);
    hold on
    plot(t_opt,v2,'LineWidth',1.5);
    legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    grid on
    % axis equal
    subplot(312)
    plot(t_opt,v3,'LineWidth',1.5);
    hold on
    plot(t_opt,v4,'LineWidth',1.5);
    grid on
    legend({'$v_3(t)$','$v_4(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    subplot(313)
    stairs(t_opt(1:N_FE:end),[results.u,nan*ones(2,1)]','LineWidth',1.5);
    legend({'$u_1(t)$','$u_2(t)$'},'interpreter','latex');
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
end
