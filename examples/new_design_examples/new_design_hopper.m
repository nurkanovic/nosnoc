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

%% Hopper OCP
% example inspired by https://github.com/KY-Lin22/NIPOCPEC and https://github.com/thowell/motion_planning/blob/main/models/hopper.jl

% The methods and time-freezing refomulation are detailed in https://arxiv.org/abs/2111.06759
%%
clear all;
close all;
clc;
import casadi.*

% delete old gif
delete hopper_simple.gif

%%
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
%%
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;  % number of stages in IRK methods
problem_options.use_fesd = 1;
problem_options.cross_comp_mode = 7;
problem_options.time_freezing = 1;
% problem_options.s_sot_max = 2;
problem_options.s_sot_min = 1;
problem_options.equidistant_control_grid = 1;
problem_options.pss_lift_step_functions = 1;
problem_options.stagewise_clock_constraint = 1;
problem_options.g_path_at_fe = 1; % evaluate path constraint on every integration step
problem_options.g_path_at_stg = 1; % evaluate path constraint on every stage point
problem_options.nonsmooth_switching_fun = 0;
problem_options.lift_complementarities = 1;
problem_options.a_n = 25;
problem_options.a_n = 1e2;
problem_options.relax_terminal_physical_time = ConstraintRelaxationMode.ELL_2;
problem_options.rho_terminal_physical_time = 1e4;


solver_options.N_homotopy = 9;
solver_options.sigma_0 = 1;
%solver_options.opts_casadi_nlp.ipopt.tol = 1e-6;
%solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
%solver_options.opts_casadi_nlp.ipopt.acceptable_iter = 3;
solver_options.opts_casadi_nlp.ipopt.max_iter = 5e3;
solver_options.complementarity_tol = 1e-7;
solver_options.print_level = 5;
solver_options.warm_start_duals = true;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% discretization
T = 1.5; % prediction horizon
N_stg = 20; % control intervals
N_FE = 3;  % integration steps per control interval

v_slip_bound = 0.001;
full_comp = 1;
x_goal = 0.7;
%% Hopper model
q = SX.sym('q', 4);
v = SX.sym('v', 4);
x = [q;v];
u = SX.sym('u', 3);
% state equations
mb = 1;  % body
ml = 0.1; % link
Ib = 0.25; % body
Il = 0.025; % link
mu = 1; % friction coeficient
mu = 0.45;
g = 9.81;

% inertia matrix
M = diag([mb + ml, mb + ml, Ib + Il, ml]);
% coriolis and graviry
C = [0;(mb + ml)*g;0;0];
% Control input matrix
B = [0, -sin(q(3));...
    0, cos(q(3));...
    1, 0;...
    0, 1];

% constraint normal
f_c_normal = [0; 1; q(4)*sin(q(3)); -cos(q(3))];
% constraint tangent;
f_c_tangent = [1; 0; q(4)*cos(q(3)); sin(q(3))];

% tangential and normal velocity of the contact point
v_tangent = f_c_tangent'*v;
v_normal = f_c_normal'*v;

% All forces
f_v = -C + B*u(1:2);

% Gap function
f_c = q(2) - q(4)*cos(q(3));

% The control u(3) is a slack for modellingo of nonslipping constraints.
ubu= [50; 50; 100];
lbu= [-50; -50; 0];

ubx = [x_goal+0.1; 1.5; pi; 0.50; 10; 10; 5; 5];
lbx =  [0; 0; -pi; 0.1; -10; -10; -5; -5];

x0 = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
x_mid = [(x_goal-0.1)/2+0.1; 0.8; 0; 0.1; 0; 0; 0; 0];
x_end = [x_goal; 0.5; 0; 0.5; 0; 0; 0; 0];

Q = diag([50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1]);
Q_terminal =diag([50; 50; 50; 50; 0.1; 0.1; 0.1; 0.1]);

Q = diag([50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1]);
Q_terminal = diag([500; 500; 300; 300; 0.1; 0.1; 0.1; 0.1]);

u_ref = [0; 0; 0];
R = diag([0.01; 0.01; 1e-5]);

% Avoid slipping motion
G_path = [v_tangent;f_c];
H_path = [u(3);u(3)];

%% interpolate refernece
x_ref = interp1([0 0.5 1],[x0,x_mid,x_end]',linspace(0,1,N_stg),'spline')'; %spline

%% Populate model
model = nosnoc.model.Cls();
problem_options.T = T;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = mu;
model.x0 = x0;

model.M = M;
model.f_v = f_v;
% gap functions
model.f_c = f_c;
model.J_tangent = f_c_tangent;
model.dims.n_dim_contact = 2;

% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
% constraint on drift velocity
%model.G_path = G_path;
%model.H_path = H_path;
% LSQ objective
model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};

%% Call nosnoc solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% read and plot results
x_res = ocp_solver.get('x');
u_res = ocp_solver.get('u');
q_opt = x_res(1:4,:);
v_opt = x_res(5:8,:);
t_opt = x_res(9,:);

% Swtiching functions: Gap function, normal velocity, tangentail velocity
c_eval = [];
c_fun = casadi.Function('c', {x}, {[f_c, v_normal, v_tangent]});
for ii = 1:length(x_res)
    % Note: Empty parameter argument as there are no global/time_varying params.
    c_eval = [c_eval;full(c_fun(x_res(1:end-1,ii)))];
end
c_eval = c_eval';
%%  plots
fig_num = 2;
plotHopperStatesControls(x_res,u_res,fig_num);
fig_num = 4;
% tagnetinal and normal velocity
plotHopperSwitchingFun(t_opt,c_eval,fig_num);

%%
if 0 
figure
% slip complement
u_normalized = [u_res(3,:),nan]./max(u_res(3,:));
slip_cc = u_normalized .*c_eval(1,1:N_FE:end);
plot(c_eval(1,1:N_FE:end));
hold on 
plot(u_normalized)
plot(slip_cc)
% yline(v_slip_bound,'r--')
% ylim([min(slip_cc)-0.01 max(min(slip_cc))+0.1])
xlabel('$t$','interpreter','latex');
ylabel('$v_{\mathrm{t}} \cdot (1-\mathrm{step}(f_c(q)))$ -tangential velocitiy','interpreter','latex');
end
%% animation
filename = 'hopper_simple.gif';

x_min = [x0(1)-0.4];
x_max = [x_ref(1,end)+0.4];
y_min = [-0.15];
y_max = [1];

x_head = q_opt(1,:);
y_head = q_opt(2,:);
x_foot = q_opt(1,:) - q_opt(4,:).*sin(q_opt(3,:));
y_foot = q_opt(2,:) - q_opt(4,:).*cos(q_opt(3,:));


x_head_ref = x_end(1);
y_head_ref = x_end(2);
x_foot_ref = x_end(1)- x_end(4).*sin(x_end(3));
y_foot_ref = x_end(2) - x_end(4).*cos(x_end(3));
q_ref = x_end(1:4);

figure(77)

for ii = 1:length(x_head)
    % The reference
    plotHopperConfiguration(q_ref,0.1);
    hold on
    plot(x_ref(1,:),x_ref(2,:),'Color',[0 0 0 0.2]);
    % The Robot
    plotHopperConfiguration(q_opt(:,ii));
    % faded trajectory
    plot(x_head(1:ii),y_head(1:ii),'Color',[0 0 1 0.5]);
    hold on
    plot(x_foot(1:ii),y_foot(1:ii),'Color',[1 0 0 0.5]);
    plotHopperScene1(x_min,x_max,y_min,y_max);

    frame = getframe(77);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',problem_options.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',problem_options.h_k(1));
    end

    if ii~=length(x_head)
        clf;
    end
end

