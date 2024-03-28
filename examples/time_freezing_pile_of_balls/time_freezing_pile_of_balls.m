clear all;
clear all;
clc;
import casadi.*
close all
%% discretization
h = 0.05;
T_sim = 1;
N_sim = T_sim/h;
N_finite_elements = 3;
% animation settings
run_animation = 1;
video_speed_up = 0.25;
%% nosnoc settings
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 1;
solver_options.print_level = 3;
solver_options.N_homotopy = 5;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 0;
solver_options.mpcc_mode = MpccMode.Scholtes_ineq;
problem_options.pss_lift_step_functions = 0;
solver_options.break_simulation_if_infeasible = 0;
%% integrator settings
model = NosnocModel();
problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = N_finite_elements;
solver_options.use_previous_solution_as_initial_guess = 1;

%% model
% dimensoon
model.a_n = 100;
model.dims.n_dim_contact = 2;
n_balls = 4;
n_q = n_balls*2; % number of positions
% parameters
m = 1;
M = m*eye(n_q);
g = 9.81;
mu = 0;
e = 0;
% Differential state
q = SX.sym('q',n_q);
v = SX.sym('v',n_q);
x = [q;v];
r = 0.2*ones(n_balls,1);
%% switching functions
q_pos = {};
for ii = 0:n_balls-1
    q_pos{ii+1} = q(2*(ii+1)-1:2*(ii+1));
end
f_c = [];
for ii = 1:n_balls-1
    for jj = ii:n_balls
        if ii~=jj
            f_c = [f_c;norm(q_pos{ii}-q_pos{jj})^2-(r(ii)+r(jj))^2];
        end
    end
end
% walls
wall_up = 2;
wall_down = 0;
wall_left = -3;
wall_right = 3;

for ii = 1:n_balls
    q_temp = q_pos{ii};
    f_walls_ii = [q_temp(2)-(wall_down+r(ii))];
    f_c = [f_c;f_walls_ii];
end
%%
% all forces without constraints (free flight dynamics)
g = 9.81;
f_g = [0;-g];
f_g = repmat(f_g,n_balls,1);
f_v = zeros(n_q,1)+f_g;
q0  = [];
for ii = 1:n_balls
    q0 = [q0;0;r(ii)+(ii-1)*2*r(ii)];
    v0 = zeros(n_q,1);
end
% x
q0(end-1) = -r(end)-2;
v0(end-1) = 20;
% y
q0(end) = q0(end)-2*r(end);
x0 = [q0;v0];
%% populate model
model.e = 0;
model.q = q;
model.v = v;
model.x = x;
model.mu_f = mu;
model.M = M;
model.f_v = f_v;
model.x0 = x0;
model.f_c = f_c;
%% Call nosnoc Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%% velocity plot
q = results.x(1:n_q,:);
v = results.x(n_q+1:end-1,:);
t_phy = results.x(end,:);
tt = 0:h:T_sim;
%% animation
hh1 = linspace(wall_left-3,wall_right+3,5);
hh2 = linspace(wall_down,wall_up,5);
q_max = max(q(:));
ind = diff(t_phy)/h>0.1;
q = results.x(1:n_q,ind);
time = t_phy(ind);
%%
delete impacting_balls
filename = 'impacting_balls';
figure('Renderer', 'painters', 'Position', [10 10 900 600])
pos = get(gcf, 'Position');
% for 4k
width = pos(3)*1.5;
height = pos(4)*1.5;
% normal
width = pos(3);
height = pos(4);
frames_gif = zeros(height, width, 1, N_sim, 'uint8');
frames_video = { };
N_frames = length(q);
for ii = 1:N_frames
    for jj = 0:n_balls-1
        plot_circle(q(2*(jj+1)-1,ii),q(2*(jj+1),ii),r(jj+1));
        hold on
    end
    axis equal
    hold on
    plot(hh1,hh1*0+wall_down,'k','linewidth',1.5);
    hold on
    text(-2,4.5,['Time: ' num2str(time(ii),'%.2f') ' s'],'interpreter','latex','fontsize',15)
    xlim([wall_left-0.5 wall_right+0.5])
    ylim([wall_down-0.5 6])
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    frame = gcf;
    set(frame,'Position',[15 15 width height])
    f = getframe(frame);
    frames_video{ii} = f;
    if ii == 1
        [frames_gif(:,:,1,ii), cmap] = imresize(rgb2ind(f.cdata, 256, 'nodither'), [height, width]);
    else
        frames_gif(:,:,1,ii) = imresize(rgb2ind(f.cdata, cmap, 'nodither'), [height, width]);
    end
    if ii~=N_frames
        clf;
    end
    if ii~=length(tt)
        clf;
    end
end
%%
% video file
pathVideoAVI = [ filename]; % filename, used later to generate mp4
% writerObj = VideoWriter(pathVideoAVI,'Uncompressed AVI');
writerObj = VideoWriter(pathVideoAVI,'MPEG-4');
target_duration = (T_sim/video_speed_up);
writerObj.FrameRate = round(N_frames/target_duration);
open(writerObj);
for ii = 1:N_frames
    writeVideo(writerObj,frames_video{ii});
end
close(writerObj); % Close the movie file
