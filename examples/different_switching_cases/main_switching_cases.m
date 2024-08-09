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
import casadi.*

%% Info
% The goal of this example is to illustrate the 4 different cases of switching that can occur in a PSS, on some simple examples.
% by chosing example_num from 1 to 4 below one can see a case of
% 1) crossing a disconituity
% 2) sliding mode
% 3) sliding on a surfce of disconinuity where a spontaneous switch can happen (nonuqnie solutions)
% 4) unique leaving of a sliding mode
switching_case = 'spontaneous_switch';
%  Options: 'crossing' 'sliding_mode', 'spontaneous_switch' , 'leave_sliding_mode', 
%% NOSNOC settings
problem_options = nosnoc.Options();

model = nosnoc.model.Pss();

% discretization parameters
% TODO: these should be problem_options
problem_options.N_sim = 10;
problem_options.T_sim = 1.5;

problem_options.N_finite_elements = 2;
problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_representation= 'differential';
problem_options.cross_comp_mode = 7;

switch switching_case
    case 'crossing'
        %% Crossing a discontinuity
        model.x0 = -1;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [2]; f_2 = [0.2];
        model.F = [f_1 f_2];
    case 'sliding_mode'
        %% Sliding mode
        model.x0 = -0.5;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [1]; f_2 = [-1];
        model.F = [f_1 f_2];
    case 'spontaneous_switch'
        %% spontaneous switch
        model.x0 = 0;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [-1]; f_2 = [1];
        model.F = [f_1 f_2];
        % implicit methods more accurate, explicit Euler enables "random"
        % leaving
        problem_options.rk_scheme = RKSchemes.EXPLICIT_RK;
        problem_options.n_s = 1;
        problem_options.N_finite_elements = 3; % set 4, 5 for different outcomes
    case 'leave_sliding_mode'
        %% leaving sliding mode in a unique way
        model.x0 = [0;0];
        x = SX.sym('x',1);
        t = SX.sym('t',1);
        model.x = [x;t];
        model.c = x;
        model.S = [-1; 1];
        f_1 = [1+t;1]; f_2 = [-1+t;1];
        model.F = [f_1 f_2];
    otherwise
        error('pick a value for example_num between 1 and 4.')
end

solver_options = nosnoc.solver.Options();
solver_options.print_level = 3;
solver_options.store_integrator_step_results = 1;

integrator = nosnoc.Integrator(model, problem_options, solver_options);

if strcmp('spontaneous_switch', switching_case)
    % initialize solver
    theta_init = [1; 0];
    integrator.set('theta', 'init', {1,1:2,1:problem_options.n_s}, theta_init);
end
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

figure
latexify_plot()
plot(t_grid,x_res(1,:))
grid on
xlabel('$t$')
ylabel('$x(t)$')
grid on

