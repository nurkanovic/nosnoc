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
clear all
clc
close all
import casadi.*

% example settings
illustrate_regions  = 1;
terminal_constraint = 1;
linear_control = 1;

%% NOSNOC settings and model
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
model = NosnocModel();
%%
problem_options.n_s = 3;
N_finite_elements = 3;

problem_options.irk_representation = 'integral';
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
solver_options.psi_fun_type = CFunctionType.SCHOLTES;
problem_options.cross_comp_mode = 7;
%problem_options.lift_complementarities = 1;

solver_options.print_level = 3;
problem_options.use_fesd = 1;
solver_options.comp_tol = 1e-9;
solver_options.sigma_N = 1e-9;
solver_options.N_homotopy = 15;
problem_options.equidistant_control_grid = 1;

problem_options.step_equilibration = 'heuristic_mean';  % heuristic_diff, heuristic_mean, l2_relaxed, l2_relaxed_scaled, direct, direct_homotopy, off
problem_options.rho_h = 1e2;

%% model equations

% Variable defintion
x1 = SX.sym('x1');
x2 = SX.sym('x2');

v1 = SX.sym('v1');
v2 = SX.sym('v2');

% Control
u1 = SX.sym('u1');
u2 = SX.sym('u2');
model.u = [u1;u2];

if linear_control
    v0  = [0;0];
    x = [x1;x2;v1;v2];
    u_max = 10;

    % dynamics
    f_11 = [-1+v1;0;u1;u2];
    f_12 = [1+v1;0;u1;u2];
    f_21 = [0;-1+v2;u1;u2];
    f_22 = [0;1+v2;u1;u2];

    % Objective
    model.f_q = 1*(v1^2+v2^2)+0*(u1^2+u2^2);
else
    u_max = 2;
    v0 = [];
    x = [x1;x2];

    % dynamics
    f_11 = [-1+u1;0];
    f_12 = [1+u1;0];
    f_21 = [0;-1+u2];
    f_22 = [0;1+u2];

    % Objective
    model.f_q = u1^2+u2^2;
end
model.x0 = [2*pi/3;pi/3;v0];
model.x = x;
problem_options.T = 4;

problem_options.N_stages = 20;
problem_options.N_finite_elements = N_finite_elements;

% Switching Functions
p = 2; a = 0.15; a1 = 0;
b = -0.05; q = 3;

c1 = x1+a*(x2-a1)^p;
c2 = x2+b*x1^q;
model.c = {c1,c2};

S1 = [1;-1];
S2 = [1;-1];
model.S = {S1,S2};

%% Modes of the ODEs layers (for all  i = 1,...,n_sys);
F1 = [f_11 f_12];
F2 = [f_21 f_22];
model.F = {F1,F2};

% constraints
model.lbu  = -u_max*ones(2,1);
model.ubu  = u_max*ones(2,1);

x_target = [-pi/6;-pi/4];
if terminal_constraint
    model.g_terminal = [x(1:2)-x_target(1:2)];
else
    model.f_q_T = 100*(x(1:2)-x_target(1:2))'*(x(1:2)-x_target(1:2));
end

%% Solve and plot
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();

u_opt = results.u;
f_opt = full(results.f);

t_grid_optimizer = [results.t_grid];
x_res_optimizer = [results.x];
%%
figure
stairs(results.t_grid,[results.h,nan])
xlabel('$t$','Interpreter','latex');
ylabel('$h_{ki}$','Interpreter','latex');
%%
fprintf('Objective value %2.4f \n',f_opt);
f_star = 6.616653254750982;
fprintf('Error %2.4e \n',norm(f_opt - f_star));


x_res_integrator = [];
t_grid_integrator = [];
t_end = 0;


tspan = [0 4/6];
y0 = [2*pi/3;pi/3;0;0];
sigma_int = solver_options.comp_tol;
tol = sigma_int/10;
options = odeset('RelTol', tol, 'AbsTol', tol/10);
x_res_integrator = [];
t_grid_integrator = [];
for ii = 1:6
    [t,y_res] = ode15s(@(t,y) ...
        ((0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([-1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+...
        ((0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;-1+y(4);u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;1+y(4);u_opt(1,ii);u_opt(2,ii)])...
        ,tspan, y0,options);

    y0 = y_res(end,:)';
    x_res_integrator = [x_res_integrator,y_res'];
    t_grid_integrator = [t_grid_integrator,t'+(ii-1)*4/6];
end
x_end = x_res_integrator(1:2,end);
fprintf('Error %2.4e \n',norm(x_target - x_end));

%%
if 1
    figure
    subplot(211)
    plot(t_grid_optimizer,x_res_optimizer(1:2,:))
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$x(t)$ - optimizer','interpreter','latex');
    subplot(212)
    plot(t_grid_integrator,x_res_integrator(1:2,:))
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$x(t)$ - integrator','interpreter','latex');
    %
    figure
    subplot(211)
    plot(t_grid_optimizer,x_res_optimizer(3:4,:))
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$ - optimizer','interpreter','latex');
    subplot(212)
    plot(t_grid_integrator,x_res_integrator(3:4,:))
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$ - integrator','interpreter','latex');

    %
    figure
    subplot(121)
    plot(x_res_optimizer(1,:),x_res_optimizer(2,:),'LineWidth',2)
    grid on
    hold on
    plot(x_target(1),x_target(2),'rx')
    if illustrate_regions
        hold on
        t2 = -5:0.01:5;
        plot(-a*(t2-a1).^p,t2,'k')
        hold on
        t1 = t2;
        plot(t1,-b*t1.^q,'k')
        grid on
        axis equal
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
    end
    subplot(122)
    plot(x_res_integrator(1,:),x_res_integrator(2,:),'LineWidth',2)
    hold on
    plot(x_target(1),x_target(2),'rx')
    grid on
    axis equal
    if illustrate_regions
        hold on
        %     p = 2; a = -0.1; a1 = 0;
        %     b = -0.05; q = 3;
        t2 = -5:0.01:5;
        plot(-a*(t2-a1).^p,t2,'k')
        hold on
        t1 = t2;
        plot(t1,-b*t1.^q,'k')
        grid on
        axis equal
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
    end
end
%%
sliding_mode_plot_for_paper
