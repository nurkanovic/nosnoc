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

%% Settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 1;
%% Time settings
problem_options.time_freezing = 1;
problem_options.time_freezing_hysteresis = 1;
problem_options.time_optimal_problem = 1;
% Time-freezing scaling / speed of time
problem_options.s_sot_max = 100;
problem_options.s_sot_min = 0.1;
problem_options.rho_sot = 1e-1;
problem_options.use_speed_of_time_variables = 1; 
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 1;
problem_options.relax_terminal_constraint = 0;
% solver settings
solver_options.comp_tol = 1e-6;
problem_options.cross_comp_mode = 4;
%solver_options.homotopy_update_slope = 0.9;
%solver_options.homotopy_update_rule = 'superlinear';
%solver_options.psi_fun_type = CFunctionType.STEFFENSEN_ULBRICH;
%solver_options.elasticity_mode = ElasticityMode.ELL_INF;
problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e4;
solver_options.opts_casadi_nlp.ipopt.tol = 1e-7;
solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-5;
solver_options.opts_casadi_nlp.ipopt.acceptable_iter = 3;
solver_options.sigma_0 = 10;
solver_options.N_homotopy = 7;
solver_options.print_level = 3;

%% Model Settings
problem_options.N_finite_elements = 3;
problem_options.N_stages = 10;
%% solve OCP
model = car_hystheresis_model_voronoi();
problem_options.T = 1;
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%% Read and plot Result
plot_results_car_hysteresis(results,problem_options,model,stats)   
%%
T_opt = results.T;
N_stages = problem_options.N_stages;
u_opt = results.u;
[tout,yout,error] = car_hysteresis_sim(u_opt,T_opt,N_stages);
