%
%    This file is part of NOSNOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
function [settings] = refine_settings_integrator(settings);
% This functions addapts the default/user provided settings so that they
% make sense for the integrator.
%% Unfold user structure
unfold_struct(settings,'caller')
clear settings;
%% Number of stages and times.
if equidistant_control_grid == 0
    couple_across_stages = 1;
end


if (time_freezing || time_optimal_problem) == 1
    time_rescaling = 1;
else
    if time_rescaling
        warning('Time rescaling makes only sense if either time freezing or time optimal problem is on. Setting time_rescaling =0.')
        time_rescaling = 0;
    else
        time_rescaling = 0;
    end
end

if time_rescaling == 0
    use_speed_of_time_variables  = 0;
end

if use_speed_of_time_variables == 0
    local_speed_of_time_variable = 0;
end

if exist('cross_complementarity_mode')
    cross_comp_mode = cross_complementarity_mode;
end

%% Impose time
% this gets active only if time freezing is activated and the user had provided impose_terminal_phyisical_time = 0.
if (time_freezing && ~impose_terminal_phyisical_time) 
    time_rescaling = 0;
    % But what if I want to get time_optimal_problem and minimize the
    % numericla time? The implementation of this should be avoided.
end
if impose_terminal_phyisical_time == 0
    warning ('impose_terminal_phyisical_time = 0 is not recommended. It means T \neq T_phy (or T_final \neq T_phy). It is only supported for nonequdistant control grids \n')
end

if time_freezing
    use_speed_of_time_variables = 0;
    local_speed_of_time_variable= 0;
    stagewise_clock_constraint = 0;
    time_freezing = 0;
    time_rescaling = 0;
end

% lifting does not make sense in integral mode
if isequal(irk_representation,'integral') 
    lift_irk_differential = 0;
end

%% Save data for output into struct
% settings = [];
names = who;
for ii = 1:length(names)
    eval([ 'settings.' names{ii} '=' names{ii} ';'])
end
end