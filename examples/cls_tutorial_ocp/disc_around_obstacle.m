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

%% Manipulation of two discs 

clear all;
close all;
clc;
import casadi.*

filename = 'discs_switch_position_obstacle.gif';
%%
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%%
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 2;  % number of stages in IRK methods

solver_options.mpcc_mode = 'elastic_ineq'; % \ell_inifnity penalization of the complementariy constraints
problem_options.use_fesd = 1;
solver_options.N_homotopy = 7;
problem_options.dcs_mode = 'CLS';
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
problem_options.g_path_at_fe = 1;
% 
solver_options.mpcc_mode = 'Scholtes_ineq';
solver_options.homotopy_update_slope = 0.5;
%solver_options.homotopy_update_rule = 'superlinear';
solver_options.N_homotopy = 100;
problem_options.g_path_at_fe = 1;
problem_options.cross_comp_mode = 7;
solver_options.sigma_0 = 1e1;
solver_options.comp_tol = 1e-6;
solver_options.opts_casadi_nlp.ipopt.max_iter = 2e3;
problem_options.gamma_h = 0.995;
%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
%solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretizatioon
N_stg = 25; % control intervals
N_FE = 2;  % integration steps per control interval
T = 4;

%% model parameters
m1 = 2;
m2 = 1;
r1 = 0.3;
r2 = 0.2;

ubx = [10; 10;10; 10; 5; 5; 5; 5];
lbx = -ubx;

q10 = [-2; 0];
q20 = [2;0];
q_target1 = q20;
q_target2 = q10;
v10 = [0;0];
v20 = [0;0];
x0 = [q10;q20;v10;v20];

u_ref = [0;0];
x_ref = [q_target1;q_target2;zeros(4,1)];

Q = diag([10;10;10;10;0*ones(4,1)]);
R = diag([0.1 0.1]);
Q_terminal = 100*Q;

ubu = [30;30];
lbu = -ubu;

%% Symbolic variables and bounds
q = SX.sym('q',4);
v = SX.sym('v',4);
u = SX.sym('u',2);

q1 = q(1:2);
q2 = q(3:4);
v1 = v(1:2);
v2 = v(3:4);

x = [q;v];
problem_options.T = T;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 1*0;
model.mu_f = 0.0;
model.a_n = 10;
model.x0 = x0;
model.dims.n_dim_contact = 2;

cv = 2;
eps = 1e-1;
f_drag = cv*[v1/norm(v1+eps);v2/norm(v2+eps)];

model.M = diag([m1;m1;m2;m2]); % inertia/mass matrix;
model.f_v = [u;...
             zeros(2,1)]-f_drag;

%% gap functions
model.f_c = [norm(q1-q2)^2-(r1+r2)^2];

%% obstacle
r_ob = 1;
q_ob = [0;0];
g_path = -[(q1(1)-q_ob(1))^2+(q1(2)-q_ob(2))^2-(r_ob+r1)^2;...
    (q2(1)-q_ob(1))^2+(q2(2)-q_ob(2))^2-(r_ob+r2)^2];
model.g_path = g_path;

%% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
%%  Objective
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
%% Call nosnoc solver
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%% read and plot results
unfold_struct(results,'base');
p1 = x(1,:);
p2 = x(2,:);
p3 = x(3,:);
p4 = x(4,:);
v1 = x(5,:);
v2 = x(6,:);
v3 = x(7,:);
v4 = x(8,:);
u_opt = u;

p1 = x_with_impulse(1,:);
p2 = x_with_impulse(2,:);
p3 = x_with_impulse(3,:);
p4 = x_with_impulse(4,:);
v1 = x_with_impulse(5,:);
v2 = x_with_impulse(6,:);
v3 = x_with_impulse(7,:);
v4 = x_with_impulse(8,:);
t_grid = results.t_with_impulse;


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
    % obstacle
    plot(r_ob*x_t+q_ob(1),r_ob*y_t+q_ob(2),'k-','LineWidth',1.5);
    plot(q_ob(1),q_ob(2),'Color',0.5*ones(3,1),'Marker','.','MarkerSize',500)
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

%%  several frames
% figure('Renderer', 'painters', 'Position', [100 100 1000 400])
figure('Renderer', 'painters', 'Position', [100 100 800 300])

x_min = min([p1,p2,p3,p4])-1;
x_max = max([p1,p2,p3,p4])+1;

tt = linspace(0,2*pi,100);
x_t = cos(tt);
y_t = sin(tt);
N_total = length(p1);
N_shots = 12;
N_skip = round(N_total/N_shots);
for jj= 1:N_shots
    if jj ~=N_shots
    ii = (jj-1)*(N_skip)+1;
    else
        ii = N_total;
    end
    subplot(2,N_shots/2,jj)
    plot(r1*x_t+p1(ii),r1*y_t+p2(ii),'b-','LineWidth',2);
    hold on
    plot(r2*x_t+p3(ii),r2*y_t+p4(ii),'r-','LineWidth',2);
    text(-1.5,2,['$t = ' num2str(round(t_grid(ii),2)) '\ s$'],'interpreter','latex');
    % obstacle
    plot(r_ob*x_t+q_ob(1),r_ob*y_t+q_ob(2),'k-','LineWidth',1.5);
    axis equal
    xlim([x_min x_max])
    ylim([x_min x_max])
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('$y$ [m]','Interpreter','latex');
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['manipulation_frames'] ])
%%
figure('Renderer', 'painters', 'Position', [100 100 800 450])
subplot(subplot(321))
plot(t_grid,p1,'LineWidth',1.5);
hold on
plot(t_grid,p2,'LineWidth',1.5);
xlim([0 T]);
legend({'$q_{1,1}(t)$','$q_{1,2}(t)$'},'interpreter','latex','location','best');
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
grid on
% axis equal
subplot(322)
plot(t_grid,p3,'LineWidth',1.5);
hold on
plot(t_grid,p4,'LineWidth',1.5);
grid on
xlim([0 T]);
legend({'$q_{2,1}(t)$','$q_{2,2}(t)$'},'interpreter','latex','location','southwest');
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');


subplot(323)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
xlim([0 T]);
legend({'$v_{1,1}(t)$','$v_{1,2}(t)$'},'interpreter','latex','location','north');
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
grid on
% axis equal
subplot(324)
plot(t_grid,v3,'LineWidth',1.5);
hold on
plot(t_grid,v4,'LineWidth',1.5);
grid on
xlim([0 T]);
legend({'$v_{2,1}(t)$','$v_{2,2}(t)$'},'interpreter','latex','location','northwest');
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
subplot(3,2,5)
stairs(t_grid_u,[u_opt,nan*ones(2,1)]','LineWidth',1.5);
legend({'$u_{1,1}(t)$','$u_{1,2}(t)$'},'interpreter','latex','location','best');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
xlim([0 T]);
set(gcf,'Units','inches');


subplot(3,2,6)
stem(results.t_grid(1:problem_options.n_s-1:end),[nan,Lambda_normal]','LineWidth',1.5);
grid on
hold on
% plot(results.t_grid(1:problem_options.n_s-1:end),[nan,lambda_normal]','LineWidth',1.5);
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}(t)$','interpreter','latex');
xlim([0 T]);
set(gcf,'Units','inches');


screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['manipulation_states'] ])

