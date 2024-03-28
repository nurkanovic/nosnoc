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

clear all
clc
close all
import casadi.*
%%
switching_case = 'sliding_mode';
%% NOSNOC settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 1;
settings.homotopy_update_slope = 0.1;
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.irk_representation= 'differential';
problem_options.irk_representation= 'integral';
problem_options.dcs_mode = 'Step';

% discretization parameters
N_sim = 1;
T_sim = 0.75;

model = NosnocModel();
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;

model.x0 = -0.50;
x = SX.sym('x',1);
model.x = x;
model.c = x;
model.S = [-1; 1];
f_1 = [1]; f_2 = [-1];
model.F = [f_1 f_2];

integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%
figure
plot(results.t_grid,results.x)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on



