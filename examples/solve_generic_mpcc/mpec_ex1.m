% Example of a small MPCC example, solved via the nonsnoc MPCC solver mpccsol. 
% The origin is C stationary, (1,0) is S-stationary
close all;
clear;
import casadi.*;
import nosnoc.solver.mpccsol;

opts = nosnoc.solver.Options();
mpccsol_opts = opts;  % TODO: remove relaxation layer
% mpccsol_opts.solver_type = '??'  % to be adressed in current PR.
mpccsol_opts.calculate_stationarity_type = true; % does nothing


x1 = SX.sym('x1');
x2 = SX.sym('x2');
% parameters
p = SX.sym('p');
x = [x1;x2];
f = (x1-1)^2+x2^3+x2^2+p;
x0 = [0;5];
p0 = 0;
lbx = [0;0];
ubx = [inf;inf];

mpcc_struct.x = x;
mpcc_struct.g = [];
mpcc_struct.p = p;
mpcc_struct.G = x1;
mpcc_struct.H = x2;
mpcc_struct.f = f;

solver = mpccsol('generic_mpcc', 'scholtes_ineq', mpcc_struct, mpccsol_opts);

mpcc_results = solver('x0', x0,...
    'lbx', lbx,...
    'ubx', ubx,...
    'p',p0);

disp(mpcc_results.x)