clear all
close all
import casadi.*
%import vdx.*

T = 5.0;
N_sim = 100;
t_step = T/N_sim;
%% Define (uncontrolled for now) projected system
model = nosnoc.pds.Model;
opts = nosnoc.Options;
x = SX.sym('x', 2);
model.x = x;
model.lbx = [-inf;-inf];
model.ubx = [inf;inf];
model.x0 = [0;0.5];
model.u = [];
model.lbu = [];
model.ubu = [];
model.u0 = [];
model.c = [x(2)+0.25;-x(2) - (x(1)+0.5)^2 + 1];
model.f_x = [x(2); -x(1)];
model.f_q = 0;
model.f_q_T = 0;

opts.T = t_step;
opts.N_stages = 1;
opts.N_finite_elements = 3;
opts.n_s = 2;
opts.cross_comp_mode = CrossCompMode.FE_FE;

prob = nosnoc.pds.Problem(model, opts);
prob.create_variables();
prob.forward_sim_constraints();
prob.generate_complementarities();
