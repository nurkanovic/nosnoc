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
classdef MpccMode
% Brief MPCC Wiki
% There are several possible MPCC Solution strategies avilable, by setting mpcc_mode to :
% 'direct' - treat complementarity conditions directly in the NLP, the bilinear term is tread as an inequality constraint.
% 'Scholtes_eq' - Smooth the complementarity conditions, Scholtes' smoothing.
% 'Scholtes_ineq' - Relax the complementarity conditions, Scholtes' relaxation.
% 'ell_1_penalty' - \ell_1 penalty, penalize the sum of all bilinear terms in the objective
% 'elastic_ineq' - \ell_infty elastic mode, upper bound all bilinear term with a positive slack, and penalize the slack in the objective.
% 'elastic_eq' - \ell_infty elastic mode, equate all bilinear term to a positive slack, and penalize the slack in the objective.
% 'elastic_two_sided' - \ell_infty, same as 'elastic_ineq' but two sided.
% 'elastic_ell_1_ineq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_eq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_two_sided' - \ell_1, elastic mode but penalize ell_1 norm of complementarities

    enumeration
        direct
        Scholtes_ineq
        Scholtes_eq
        ell_1_penalty
        elastic_ineq
        elastic_eq
        elastic_two_sided
        elastic_ell_1_ineq
        elastic_ell_1_eq
        elastic_ell_1_two_sided
        lifting
        custom
    end

    properties(Constant)
        elastic_ell_1 = ["elastic_ell_1_ineq"
                         "elastic_ell_1_eq"
                         "elastic_ell_1_two_sided"];
        elastic = ["elastic_ineq"
                   "elastic_eq"
                   "elastic_two_sided"];
    end
    
end

