% Example of a small MPCC example, solved via the nonsnoc MPCC solver mpccsol. 
% The origin is C stationary, (1,0) is S-stationary
close all;
clear;
import casadi.*;
import nosnoc.solver.*;

mpccsol_opts = nosnoc.solver.Options();  


x1 = SX.sym('x1');
x2 = SX.sym('x2');
% parameters
p = SX.sym('p');
x = [x1;x2];
f = (x1-1)^2+x2^3+x2^2+p;
x0 = [0;0];
p0 = 0;
lbx = [0;0];
ubx = [inf;inf];

mpcc = vdx.problems.Mpcc();
mpcc.w.x = {x, lbx, ubx, x0};
mpcc.p.p = {p, p0};
mpcc.G.G = {x1};
mpcc.H.H = {x2};
mpcc.f = f;

mpcc_method1 = 'scholtes_ineq';
mpcc_method2 = 'kanzow_schwartz_ineq';
mpcc_method3 =  nosnoc.solver.MpccMethod.SCHOLTES_EQ;

mpccsol_opts.homotopy_steering_strategy = "DIRECT";
% the parameter relaxation/smoothing parameter s is minmized in the objective, its penality (1/sigma) is steered outside, cf. Table 1 in https://arxiv.org/pdf/2312.11022.pdf 
% mpccsol_opts.homotopy_steering_strategy = "DIRECT"; 
% the parameter relaxation/smoothing parameter sigma is steered outside the optimization , cf. Table 1 in https://arxiv.org/pdf/2312.11022.pdf

mpccsol_opts.calculate_stationarity_type  = 1;

%% create solver
solver = mpccsol('nosnoc_solver', mpcc_method1, mpcc, mpccsol_opts);
solver.generate_c_solver('');
