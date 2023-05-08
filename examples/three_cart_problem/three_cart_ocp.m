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
% example taken from https://github.com/KY-Lin22/NIPOCPEC

%% 
clear all;
close all;
clc;
import casadi.*

% delete old gif
delete three_carts.gif

%%
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 1;  % number of stages in IRK methods

settings.mpcc_mode = 'elastic_ineq'; % \ell_inifnity penalization of the complementariy constraints
settings.N_homotopy = 7;
settings.opts_casadi_nlp.ipopt.max_iter = 1e3;
settings.print_level = 3;
settings.comp_tol = 1e-9;
settings.time_freezing = 1;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretizatioon
N_stg = 20; % control intervals
N_FE = 3;  % integration steps per control intevral

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;
k1 = 1*0;
k2 = 1*0;
k3 = 1*0;
c_damping = 2;

ubx = [10; 15; 15; 5; 5; 5]; 
lbx = [-15; -15; -10; -5; -5; -5];            
        
x0 = [-3; 0; 3; 0; 0; 0];
x_ref = [-5; 0; 5; 0; 0; 0];
u_ref = [0; 0; 0];

Q = diag([10; 1; 10; 0.1; 0.1; 0.1]);
Q_terminal = 100*Q;
R = diag([0.1; 0.1; 0.1]);

u_max = 20;
u_min = -20;

lbu = [0*u_min; u_min; 0*u_min];
ubu = [0*u_max; u_max; 0*u_max];


%% Symbolic variables and bounds
q = SX.sym('q',3); v = SX.sym('v',3); 
u = SX.sym('u',3);
x = [q;v];
model.T = 4;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.0;
model.a_n = 20;
model.x0 = x0; 

% model.f = [1/m1*(u(1)-c_damping*v(1)-k1*q(1));...
%            1/m2*(u(2)-c_damping*v(2)+k2*(q(1)-q(2)));...
%            1/m3*(u(3)-c_damping*v(3)+k3*q(2))];
% dot{v} r.h.s.

model.M = diag([m1;m2;m3]); % inertia/mass matrix;
model.f_v = [(u(1)-c_damping*v(1)-k1*q(1));...
           (u(2)-c_damping*v(2)+k2*(q(1)-q(2)));...
           (u(3)-c_damping*v(3)+k3*q(2))];

% gap functions
model.f_c = [q(2) - q(1) - 0.5*cart_width2 - 0.5*cart_width1;...
           q(3) - q(2) - 0.5*cart_width3 - 0.5*cart_width2];
model.n_dim_contact = 2;


% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
% terminal cost
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
% model.g_terminal = [x-x_ref];
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
% [results] = polish_homotopy_solution(model,settings,results,stats.sigma_k); % (experimental, projects and fixes active set at solution and solves NLP)
%% read and plot results
unfold_struct(results,'base');
p1 = x_opt(1,:);
p2 = x_opt(2,:);
p3 = x_opt(3,:);
v1 = x_opt(4,:);
v2 = x_opt(5,:);
v3 = x_opt(6,:);
t_opt = x_opt(7,:);

u1_opt = u_opt(1,:);
u2_opt = u_opt(2,:);
u3_opt = u_opt(3,:);


%% animation
% figure('Renderer', 'painters', 'Position', [100 100 1000 400])
figure(1)
filename = 'three_carts.gif';
carts_appart = 2;
x_min =min([p1,p2,p3])-2.5;
x_max = max([p1,p2,p3])+2.5;
cart_height = 2;

carts_appart = 1.5*1;
for ii = 1:length(p1)
    % cart 1
    xp = [p1(ii)-cart_width1/2 p1(ii)+cart_height/2 p1(ii)+cart_height/2 p1(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.8)
    hold on
    % cart 2
    xp = [p2(ii)-cart_width2/2 p2(ii)+cart_height/2 p2(ii)+cart_height/2 p2(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.8)

    % cart 3
    xp = [p3(ii)-cart_width3/2 p3(ii)+cart_height/2 p3(ii)+cart_height/2 p3(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.8)

    if k2 >0
    spring1 = linspace(min(p1(ii),p2(ii)),max(p1(ii),p2(ii)),5);
    plot(spring1,spring1*0+cart_height/2,'ko-','LineWidth',1);
    end

    % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.05,'EdgeColor','none')
    hold on
    % cart 2
    xp = [x_ref(2)-cart_width2/2 x_ref(2)+cart_height/2 x_ref(2)+cart_height/2 x_ref(2)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.05,'EdgeColor','none')

    % cart 3
    xp = [x_ref(3)-cart_width3/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.05,'EdgeColor','none')


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

    if ii~=length(p1)
        clf;
    end
end

%%
figure('Renderer', 'painters', 'Position', [100 100 1400 600])
subplot(131)
plot(t_opt,p1,'LineWidth',1.5);
hold on
plot(t_opt,p2,'LineWidth',1.5);
plot(t_opt,p3,'LineWidth',1.5);
% axis equal
grid on
legend({'$p_1(t)$','$p_2(t)$','$p_3(t)$'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$p$','interpreter','latex');
% axis equal
subplot(132)
plot(t_opt,v1,'LineWidth',1.5);
hold on
plot(t_opt,v2,'LineWidth',1.5);
plot(t_opt,v3,'LineWidth',1.5);
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

subplot(133)
stairs(t_opt(1:N_FE:end),[u1_opt,nan],'LineWidth',1.5);
hold on
stairs(t_opt(1:N_FE:end),[u2_opt,nan],'LineWidth',1.5);
stairs(t_opt(1:N_FE:end),[u3_opt,nan],'LineWidth',1.5);
legend({'$u_1(t)$','$u_2(t)$','$u_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u$','interpreter','latex');
