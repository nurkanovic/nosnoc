clc; close all; clear all;
import casadi.*
%% settings
% collocation settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
% problem_options.time_freezing_time_rescaling = 1;
problem_options.use_speed_of_time_variables =  1; 
problem_options.local_speed_of_time_variable = 1;  
problem_options.stagewise_clock_constraint = 1;
problem_options.time_freezing = 1;
solver_options.N_homotopy = 6;
problem_options.N_stages = 40;
problem_options.N_finite_elements = 4;
% solver_options.homotopy_update_rule = 'superlinear';
solver_options.print_level = 5;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
%% Generate Model
% angulary velocity of reference
omega = -2*pi; % no impacts
omega = -3*pi; % impacts
model = elastic_ball_in_box_model(omega);
%% Solve OCP via nosnoc
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%% Read and plot Result
unfold_struct(results,'base');
unfold_struct(model,'base');
plot_results_elastic_ball_in_box
