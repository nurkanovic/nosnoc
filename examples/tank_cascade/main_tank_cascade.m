import casadi.*
clear all
clc 
close all
%% Descrpition
% this is an optimal control example from 
% Baumrucker, Brian T., and Lorenz T. Biegler. "MPEC strategies for optimization of a class of hybrid dynamic systems." Journal of Process Control 19.8 (2009): 1248-1256.

%% settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;   
problem_options.n_s = 1;                
% solver_options.N_homotopy = 8;
solver_options.homotopy_update_rule = 'superlinear';
%% Generate Model
model = tank_cascade();
%% Discretization parameters
problem_options.N_stages = 100; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval
problem_options.T = 100;    % yime horizon
%% Solve OCP via nosnoc
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();
%% Read and plot Result 
figure
subplot(211)
plot(results.t_grid,results.x);
hold on
xlabel('$t$','interpreter','latex');
ylabel('$L(t)$','interpreter','latex');
subplot(212)
stairs(results.t_grid_u,[nan*ones(model.dims.n_u,1),results.u]');
hold on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');














