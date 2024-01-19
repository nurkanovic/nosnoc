% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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

import casadi.*

%% Leyffer 2007, example 1, origin is C stationary, (0,1) and (1,0) are S-stationary

x1 = SX.sym('x1');
x2 = SX.sym('x2');

x = [x1;x2];

% parameters
p = SX.sym('p');
p0 = 2;

f = x1^2+x2^2-p*(x1+x2)+2;
comp1 = x1;
comp2 = x2;

x0 = [3.4;7];
lbx = [0;0];
ubx = [inf;inf];

g = [];
lbg = [];
ubg = [];

mpcc.w = x;
mpcc.w0 = x0;
mpcc.lbw = lbx;
mpcc.ubw = ubx;

mpcc.f = f;

mpcc.g = g;
mpcc.lbg = lbg;
mpcc.ubg = ubg;

mpcc.p = p;
mpcc.p0 = p0;

mpcc.G = x1;
mpcc.H = x2;

solver_options = NosnocSolverOptions();
%%
solver = NosnocSolver(mpcc, solver_options);

[results, stats] = solver.solve();
