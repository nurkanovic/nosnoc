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

% Fixed friction force
F_friction = 2;

% model
model = get_cart_pole_with_friction_model(true, F_friction);
x_ref = [0; 180/180*pi; 0; 0]; % target position

% Discretization options
problem_options = nosnoc.Options();
problem_options.T = 5;  % Time horizon
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 3;
problem_options.dcs_mode = 'Stewart';
problem_options.N_stages = 20; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control interval
problem_options.cross_comp_mode = 7;

% solver options
solver_options = nosnoc.solver.Options();
solver_options.N_homotopy = 15;
solver_options.complementarity_tol = 1e-13;
solver_options.sigma_N = 1e-13;

% other linear solvers require installation, check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions
% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

% evaluate
x_opt = ocp_solver.get("x");
distance_to_target = abs(x_ref-x_opt(:,end));
disp(['final difference to desired angle: ', num2str(distance_to_target(2), '%.3e'), ' rad'])

% visualtize
plot_cart_pole_trajectory(ocp_solver, problem_options.h_k(1), x_ref)
