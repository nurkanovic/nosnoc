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

% This file is part of NOSNOC.%
clear all; clc; close all
import casadi.*
%%
% An optimal control problem from: 
% Optimal control of systems with discontinuous differential equations 
% April 2012, Numerische Mathematik 114(4):653-695
% DOI: 10.1007/s00211-009-0262-2
% David E. Stewart Mihai Anitescu

%% Problem parameters
path_constraint = 'track';
track_width = 0.5;
omega = 1;
chicane_tightness = 1;
chicane_width = 2;

%% NOSNOC settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 4; 
solver_options.sigma_0 = 1e0;
problem_options.use_fesd = 1;
problem_options.cross_comp_mode = 1;
problem_options.T_final_max = 5*pi;
problem_options.T_final_min = 2;
%solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
problem_options.dcs_mode = 'Step';
problem_options.g_path_at_fe = 1;
problem_options.g_path_at_stg = 1;
problem_options.time_optimal_problem = 1;
solver_options.mpcc_mode = 'Scholtes_ineq';
solver_options.opts_casadi_nlp.ipopt.max_iter = 1000;
solver_options.print_level = 3;


%% model equations
model = NosnocModel();
model.x0 = zeros(5,1);
problem_options.T = 5; 
problem_options.N_stages = 50; 
problem_options.N_finite_elements = 2;

%% Variable defintion
% Declare model variables
qx = SX.sym('qx'); qy = SX.sym('qy');
vx = SX.sym('vx'); vy = SX.sym('vy');
alpha = SX.sym('alpha');
tangent = [cos(alpha);sin(alpha)];
normal = [-sin(alpha);cos(alpha)];
q = [qx;qy];
v = [vx;vy];
x = [q;v;alpha];
model.x = x;
model.lbx = [-10;-10;-inf;-inf;-inf];
model.ubx = [10;10;inf;inf;inf];

%%  control
a = SX.sym('a');
s = SX.sym('s');
u = [a;s];
u_max  = 2;
model.u = u;
model.lbu = -u_max*ones(2,1);
model.ubu = u_max*ones(2,1);
model.u0 = [u_max;0];
%% modes of the PSS
Friction_max = 4;
f_1 = [v;a*tangent+Friction_max*normal;s*(tangent'*v)];
f_2 = [v;a*tangent-Friction_max*normal;s*(tangent'*v)];
% Switching functions and modes
model.c = normal'*v;
model.S = [1;-1];
model.F = [f_1 f_2];
%%  general nonlinear constraints
switch path_constraint 
    case 'none'
        q_target = [3*pi;0];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = 0*xx;
    case 'linear'
        model.g_path = qy-(qx);
        q_target = [3*pi;3*pi];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = xx;

    case 'nonlinear'
        model.g_path = qy-sin(omega*qx);
        q_target = [3*pi;sin(omega*3*pi)];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = sin(omega*xx);
    case 'track'
        arg1 = qx-pi;
        arg2 = qx-2*pi;
        sig = 1e-1;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        model.g_path = qy-(sin(qx)*(1-step1)+(pi-qx)*step1*(1-step2)+(-pi-sin(qx))*step2);
        q_target = [3*pi;-pi];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = sin(xx).*(1-0.5.*(1+tanh((xx-pi)/sig))) +...
             (pi-xx).*0.5.*(1+tanh((xx-pi)/sig)).*(1-0.5.*(1+tanh((xx-2*pi)/sig)))+...
             (-pi-sin(xx)).*0.5.*(1+tanh((xx-2.*pi)/sig));
    case 'chicane'
        
        q_target = [10;2*chicane_width];
        model.g_path = qy-((chicane_width)+chicane_width*tanh(chicane_tightness*(qx-q_target(1)/2)));
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
end
if ~isequal(path_constraint,'none')
    model.g_path_lb = [-track_width];
    model.g_path_ub = [+track_width];
end
%% objective
if ~problem_options.time_optimal_problem 
    model.f_q = u'*u;
    model.T = 5;
else
    model.f_q = 0.01*u'*u;
end

% Terminal Constraint
model.g_terminal = [q-q_target];

% intitial guess for fun
alpha = [atan2(diff(yy),diff(xx)),0];
x_guess = vertcat(xx,yy, zeros(2,problem_options.N_stages), alpha);
%% Solve
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
%solver.set('x',mat2cell(x_guess', [ones(problem_options.N_stages,1)]));
[results,stats] = solver.solve();
fprintf('Objective values is: %2.4f \n',full(results.f));
%% plot
plot_results_ms(model, problem_options, results,...
    q_target, path_constraint, track_width, omega, chicane_tightness, chicane_width)
