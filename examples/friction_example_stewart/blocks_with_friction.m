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

function [model] = blocks_with_friction()
import casadi.*

model = nosnoc.model.Pss();
%% Initial value
model.x0 = [-1;1;-1;-1;1;1;0];
model.u0 = 0; % guess for control variables
%% Numer of ODE layers
% n_sys = 3;% number of Cartesian products in the model ("independet switches"), we call this layer
% % number of modes in every simplex
% m_1 = 2;
% m_2 = 2;
% m_3 = 2;
% n_f_sys = [m_1 m_2 m_3];

%% Variable defintion
% differential states
q1 = SX.sym('q1');
q2 = SX.sym('q2');
q3 = SX.sym('q3');
v1 = SX.sym('v1');
v2 = SX.sym('v2');
v3 = SX.sym('v3');
t = SX.sym('t');

q = [q1;q2;q3];
v = [v1;v2;v3];
model.x = [q;v;t];

%% Control
% u = SX.sym('u');
% n_u = 1;  % number of parameters,  we model it as control variables and merge them with simple equality constraints
% 
% % Guess and Bounds
% u0 = 0;
% lbu  = -20*0;
% ubu  = 20*0;

%% Switching Functions
% every constraint function corresponds to a simplex (note that the c_i might be vector valued)
c1 = v1;
c2 = v2;
c3 = v3;
% sign matrix for the modes
S1 = [1;-1];
S2 = [1;-1];
S3 = [1;-1];
% discrimnant functions
model.S = {S1,S2,S3};
model.c = {c1,c2,c3};


%% Modes of the ODEs layers (for all  i = 1,...,n_sys);
% part independet of the nonsmoothness
F_external = 0; % external force, e.g., control
F_input = 10; % variable force exicting
f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_input*cos(pi*t);...
    ]/6;

f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_external+F_input*(1*0+1*cos(pi*t));...
    1];
%
% for c1, h1,
f_11 = f_base+[0;0;0;-0.3;0;0;0];
f_12 = f_base+[0;0;0;+0.3;0;0;0];
% for c2, h2
f_21 = [0;0;0;0;-0.3;0;0];
f_22 = [0;0;0;0;0.3;0;0];
% for c3, h3
f_31 = [0;0;0;0;0;-0.3;0];
f_32 = [0;0;0;0;0;0.3;0];
% unfold_struct(model,'base');
% in matrix form
F1 = [f_11 f_12];
F2 = [f_21 f_22];
F3 = [f_31 f_32];

model.F = {F1 F2 F3};

%% Objective
% f_q = 0*u^2 + 1*v'*v;
% f_q_T = 10*v'*v;
end
