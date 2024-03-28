close all;
clear all;
import casadi.*;
import nosnoc.solver.mpccsol;

mpccsol_opts = struct;
opts = nosnoc.solver.Options();
mpccsol_opts.relaxation = opts;
mpcc = vdx.problems.Mpcc();

x = SX.sym('x', 2);
mpcc.w.x(0) = {x,0, inf};
mpcc.f = norm(x- [-1;0])^2;
mpcc.G.x(0) = {x(1)};
mpcc.H.x(0) = {x(2)};

mpcc_struct = struct;
mpcc_struct.x = mpcc.w.w;
mpcc_struct.g = mpcc.g.w;
mpcc_struct.p = mpcc.p.w;
mpcc_struct.G = mpcc.G.w;
mpcc_struct.H = mpcc.H.w;
mpcc_struct.f = mpcc.f;

solver = mpccsol('generic_mpcc', 'relaxation', mpcc_struct, mpccsol_opts);

mpcc_results = solver('x0', mpcc.w.init,...
    'lbx', mpcc.w.lb,...
    'ubx', mpcc.w.ub);

