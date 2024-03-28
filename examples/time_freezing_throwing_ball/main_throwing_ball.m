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
%% model parameters
e = 0.9; u_max = 9; beta = 0.0; 
%% NOSNOC settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.time_freezing = 1; 
problem_options.n_s = 3; 
solver_options.mpcc_mode = 'Scholtes_ineq';
solver_options.homotopy_update_rule = 'superlinear';
%solver_options.print_level=1;
%problem_options.step_equilibration = 'direct_homotopy';
q_target = [4;0.5];

model = NosnocModel();
problem_options.T = 4; 
problem_options.N_stages = 20; 
problem_options.N_finite_elements = 3;
% model equations
q = SX.sym('q',2);
v = SX.sym('v',2); 
u = SX.sym('u',2);
model.x0 = [0;0.5;0;0];
model.x = [q;v]; model.u = u; model.e = e ;
model.f_c = q(2);
model.f_v = [u-[0;9.81]-beta*v*sqrt(v(1)^2^2+v(2)^2+1e-3)]; 
% Objective and constraints
model.f_q = u'*u; model.f_q_T = 100*v'*v;
model.g_path = u'*u-u_max^2;
model.g_terminal = q-[q_target];

mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%%
plot_result_ball(model,problem_options,solver_options,results,stats)
fprintf('Objective values is: %2.4f \n',full(results.f));
fprintf('Final time is: %2.4f \n',full(results.T));
if isempty(results.T)
    results.T_opt = results.t_grid(end);
end
[tout,yout,error] = bouncing_ball_sim(results.u,results.T,problem_options.N_stages,model.x0(1:4),beta,0.9,q_target);
