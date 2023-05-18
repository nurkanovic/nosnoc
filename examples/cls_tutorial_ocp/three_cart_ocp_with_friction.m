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

%% Three cart manipulation example
clear all;
close all;
clc;
import casadi.*

%%
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
settings.n_s = 2;  % number of stages in IRK methods
settings.dcs_mode = 'CLS';
settings.opts_casadi_nlp.ipopt.max_iter = 2e3;
settings.cross_comp_mode = 1;


settings.sigma_0 = 1e1;
settings.homotopy_update_slope = 0.2;
settings.homotopy_update_rule = 'superlinear';
settings.N_homotopy = 6;
settings.friction_model = "Conic";
settings.conic_model_switch_handling = "Abs";
settings.gamma_h = 0.995;

% settings.sigma_0 = 1e1;
% settings.mpcc_mode = "elastic_ineq";
% settings.elastic_scholtes = 1;
%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretizatioon
N_stg = 20; % control intervals
N_FE = 2;  % integration steps per control intevral
T = 6;

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;

M = diag([m1, m1, m2, m2, m3, m3]);

ubx = [10;  inf; 10;    inf;  10;  inf;  10; 10; 10; 10; 10; 10]; 
lbx = [-10; -inf; -10; -inf; -10; -inf; -10; -10; -10; -10; -10;-10];            
        
x0 = [ -3; 1; 0; 1;  3; 1; ...
        0; 0; 0; 0; 0; 0];
x_ref = [-7; 0; 0; 0; 5; 0;...
          0; 0; 0; 0; 0; 0];
u_ref = 0;

Q = diag([10; 0.1; 1; 0.1; 10; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]);
Q_terminal = 200*Q;
R = 0.1;

u_max = 30;
u_min = -30;

%% Symbolic variables and bounds
g = 9.81;
q = SX.sym('q',6);
v = SX.sym('v',6); 
u = SX.sym('u',1);
x = [q;v];

q1 = q(1:2);
q2 = q(3:4);
q3 = q(5:6);

model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.x0 = x0; 

model.M = M;
model.f_v = [ 0;...
             -m1*g;
              u;...
             -m2*g;
              0;...
             -m3*g];

% gap functions
model.f_c = [q2(1) - q1(1) - 0.5*cart_width2 - 0.5*cart_width1;...
             q3(1) - q2(1) - 0.5*cart_width3 - 0.5*cart_width2;...
             q1(2)-cart_width1/2;...
             q2(2)-cart_width2/2;...
             q3(2)-cart_width3/2;...
             ];

J_tangent = [0   0 1 0 0 ;...
             -1 -1 0 0 0;...
              0  0 0 1 0;...
              1  0 0 0 0;...
              0  0 0 0 1;...
              0  1 0 0 0];

model.J_tangent = J_tangent;
model.D_tangent = [J_tangent, -J_tangent];

model.e =  [0.0 1.0 0.0 0.0 0.0];
model.mu = [0.0 0.0 0.3 0.2 0.2];

% box constraints on controls and states
model.lbu = u_min;
model.ubu = u_max;
model.lbx = lbx;
model.ubx = ubx;
% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
% model.g_terminal = [x-x_ref];
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
lambda_guess = {};
for ii = 1:N_stg
    lambda_guess{ii} = [0;0;g;g;g];
end
solver.set('lambda_normal',lambda_guess')
[results,stats] = solver.solve();
%% read and plot results
unfold_struct(results,'base');
p1x = results.x(1,:);
p2x = results.x(3,:);
p3x = results.x(5,:);
p1y = results.x(2,:);
p2y = results.x(4,:);
p3y = results.x(6,:);
v1x = results.x(7,:);
v2x = results.x(9,:);
v3x = results.x(11,:);
v1y = results.x(8,:);
v2y = results.x(10,:);
v3y = results.x(12,:);
t_opt = results.t_grid;

%% animation
% figure('Renderer', 'painters', 'Position', [100 100 1000 400])
figure(1)
filename = 'three_carts_with_friction.gif';
carts_appart = 2;
x_min = min(x_ref)-2.5;
x_max = max(x_ref)+2.5;
cart_height = 2;

carts_appart = 1.5*1;
for ii = 1:length(p1x)
    % cart 1
    xp = [p1x(ii)-cart_width1/2 p1x(ii)+cart_height/2 p1x(ii)+cart_height/2 p1x(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.8)
    hold on
    % cart 2
    xp = [p2x(ii)-cart_width2/2 p2x(ii)+cart_height/2 p2x(ii)+cart_height/2 p2x(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.8)

    % cart 3
    xp = [p3x(ii)-cart_width3/2 p3x(ii)+cart_height/2 p3x(ii)+cart_height/2 p3x(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.8)

    % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.10)
    hold on
    % cart 2
    xp = [x_ref(3)-cart_width2/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.10)

    % cart 3
    xp = [x_ref(5)-cart_width3/2 x_ref(5)+cart_height/2 x_ref(5)+cart_height/2 x_ref(5)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.10)

    % ground     
    xp = [x_min x_max x_max x_min ];
    yp = [-1 -1 0 0];
    patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');
    
    axis equal
    xlim([x_min x_max])
    ylim([-0.75 3.5])
    pause(solver.model.h_k);
   
    frame = getframe(1);
    im = frame2im(frame);

    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',solver.model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',solver.model.h_k(1));
    end

    if ii~=length(p1x)
        clf;
    end
end

%%
figure('Renderer', 'painters', 'Position', [100 100 1400 600])
subplot(131)
plot(t_opt,p1x,'LineWidth',1.5);
hold on
plot(t_opt,p2x,'LineWidth',1.5);
plot(t_opt,p3x,'LineWidth',1.5);
% axis equal
grid on
legend({'$p_1(t)$','$p_2(t)$','$p_3(t)$'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$p$','interpreter','latex');
% axis equal
subplot(132)
plot(t_opt,v1x,'LineWidth',1.5);
hold on
plot(t_opt,v2x,'LineWidth',1.5);
plot(t_opt,v3x,'LineWidth',1.5);
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

subplot(133)
stairs(t_opt(1:N_FE:end),[results.u,nan],'LineWidth',1.5);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u$','interpreter','latex');
