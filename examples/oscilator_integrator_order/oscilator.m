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

function [model] = oscilator(model_in)

import casadi.*
model = nosnoc.model.Pss();    
%% Time horizon
if ~isempty(model_in)
    unfold_struct(model_in,'caller');
else
    disp('empty struct input');
    smooth_model = 0;
end
%% Model parameters
omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end
%% Inital Value
model.x0 = [exp(-1);0];
% x0 = [2*exp(-1);0];
if ~exist('R_osc')
    R_osc = 1;
end

%% Variable defintion
% x1 = MX.sym('x1');
% x2 = MX.sym('x2');
x1 = SX.sym('x1');
x2 = SX.sym('x2');
model.x = [x1;x2];
% every constraint function corresponds to a simplex (note that the c_i might be vector valued)
c = x1^2+x2^2-R_osc^2;
% sign matrix for the modes
model.S = [1;-1];
model.c = [c];

f_11 = A1*x;
f_12 = A2*x;
% in matrix form
model.F = [f_11 f_12];

end

