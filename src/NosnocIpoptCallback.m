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

classdef NosnocIpoptCallback < casadi.Callback
    properties
        nx
        ng
        np

        nlp
        settings

        mpcc
    end
    methods
        function obj = NosnocIpoptCallback(name, mpcc, nlp, settings, nx, ng, np)
            %se...@casadi.Callback();
            obj.nx = nx;
            obj.ng = ng;
            obj.np = np;
            obj.mpcc = mpcc;
            obj.nlp = nlp;
            obj.settings = settings;

            opts.input_scheme = casadi.nlpsol_out();

            opts.output_scheme = 'ret';

            obj.construct(name,opts);

            disp('Callback instance created')

        end

        function v=get_n_in(obj)
            v=casadi.nlpsol_n_out;
        end

        function v = get_n_out(obj)
            v = 1;
        end

        function v = get_sparsity_in(obj, i)
            n = casadi.nlpsol_out(i);
            if n=='f'
                v =  casadi.Sparsity.scalar();
            elseif strcmp(n,'x') || strcmp(n,'lam_x')
                v = casadi.Sparsity.dense(obj.nx);
            elseif strcmp(n,'g') || strcmp(n,'lam_g')
                v = casadi.Sparsity.dense(obj.ng);
            elseif strcmp(n,'lam_p')
                v = casadi.Sparsity.dense(obj.np);
            else
                v = casadi.Sparsity(0,0);
            end
        end

        function [returncode] = eval(obj, arg)
            results = struct();
            results.nlp_results = struct();
            results.nlp_results.x = arg{1};
            results.nlp_results.f = arg{2};
            results.nlp_results.g = arg{3};
            obj.nlp.w.res = full(arg{1});
            obj.nlp.g.eval = full(arg{3});
            obj.nlp.f_result = full(arg{2});
            % TODO @anton add extra arguments option in the form of a cell array
            obj.settings.ipopt_callback(obj.mpcc, obj.nlp, obj.settings, [], results);
            returncode = {0};
        end

        function delete(obj)
            nlp = [];
            settings = [];
            solver = [];
        end
    end
end
