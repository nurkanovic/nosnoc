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
solver_options = NosnocSolverOptions();
model = NosnocModel();
%%
problem_options.n_s = 3;
N_finite_elements = 3;

problem_options.irk_representation = 'integral';
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.dcs_mode = 'Step';
%solver_options.psi_fun_type = CFunctionType.SCHOLTES;
%solver_options.mpcc_mode = MpccMode.elastic_ineq;
problem_options.cross_comp_mode = 7;
%problem_options.lift_complementarities = 1;

solver_options.print_level = 5;
problem_options.use_fesd = 1;
solver_options.comp_tol = 1e-5;
solver_options.sigma_N = 1e-5;
solver_options.N_homotopy = 15;
problem_options.equidistant_control_grid = 1;

problem_options.step_equilibration = 'heuristic_mean';  % heuristic_diff, heuristic_mean, l2_relaxed, l2_relaxed_scaled, direct, direct_homotopy, off
problem_options.rho_h = 1e2;

%% model equations
r = 1;
% Variable defintion
x_act = SX.sym('x_act', 2);
x_square = SX.sym('x_square', 2);

% Control
u1 = SX.sym('u1');
u2 = SX.sym('u2');
u = [u1;u2];
model.u = u;

x = [x_act;x_square];
u_max = 5;

% dynamics
f_11 = [u1;u2;0;0];
f_21 = [u1;u2;0;0];
f_31 = [u1;u2;0;0];
f_41 = [u1;u2;0;0];
f_12 = f_11 + [0;1;0;-1];
f_22 = f_21 + [1;0;-1;0];
f_32 = f_31 + [0;-1;0;1];
f_42 = f_41 + [-1;0;1;0];

% switching
c1 = x(2) - x(1) - x(4) + x(3);
c2 = x(2) + x(1) - x(3) - x(4);
c3 = x(2) - x(4) - r;
c4 = x(1) - x(3) - r;
c5 = -x(2) + x(4) - r;
c6 = -x(1) + x(3) - r;
model.c = {[c1;c2;c3;c4;c5;c6]};

S = [1,1,1,0,0,0;...
    -1,1,0,1,0,0;...
    -1,-1,0,0,1,0;...
    1,-1,0,0,0,1;...
    1,1,-1,0,0,0;...
    -1,1,0,-1,0,0;...
    -1,-1,0,0,-1,0;...
    1,-1,0,0,0,-1];
model.S = {S};


% Objective
model.f_q = 0.001*sum(u.^2);
model.f_q_T = 1000*(x(3)^2 + x(4)^2);
model.x0 = [1;3;-1;0];
model.x = x;
problem_options.T = 5;

problem_options.N_stages = 20;
problem_options.N_finite_elements = N_finite_elements;

%% Modes of the ODEs layers (for all  i = 1,...,n_sys);
F = [f_11 f_21 f_31 f_41 f_12 f_22 f_32 f_42];
model.F = {F};

% constraints
model.lbu  = -u_max*ones(2,1);
model.ubu  = u_max*ones(2,1);


%% Solve and plot
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();

u_opt = results.u;
f_opt = full(results.f);

t_grid_optimizer = [results.t_grid];
x_res_optimizer = [results.x];
