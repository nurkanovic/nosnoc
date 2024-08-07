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

function plot_cart_pole_trajectory(results, time_step, x_ref)
    link_length = 1;

    % extract results
    q1_opt = results.x(1,:);
    q2_opt = results.x(2,:);
    v1_opt = results.x(3,:);
    v2_opt = results.x(4,:);
    t_grid = results.t_grid;
    t_grid_u = results.t_grid_u;
    u_opt = results.u;

    %% Animation
    filename = 'cart_pole_with_friction.gif';
    figure('Renderer', 'painters', 'Position', [100 100 1200 600])
    cart_center = 0.125;
    cart_width1 = 0.25;
    cart_height = cart_center*2;
    pole_X = [q1_opt',q1_opt'+(link_length)*cos(q2_opt'-pi/2)];
    pole_Y = [cart_center+0*q1_opt',cart_center+link_length*sin(q2_opt'-pi/2)];
    x_min =-3;
    x_max = 3;
    for ii = 1:length(q1_opt)
        % pole
        plot(pole_X(ii,:),pole_Y(ii,:),'k','LineWidth',3);
        hold on
        % tail
        plot(pole_X(1:ii,2),pole_Y(1:ii,2),'color',[1 0 0 0.5],'LineWidth',0.5);
        % cart
        xp = [q1_opt(ii)-cart_width1/2 q1_opt(ii)+cart_height/2 q1_opt(ii)+cart_height/2 q1_opt(ii)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'k','FaceAlpha',0.8)

        % targent
        % pole
        plot([x_ref(1),x_ref(1)+(link_length)*cos(x_ref(2)-pi/2)],...
            [cart_center+0,cart_center+link_length*sin(x_ref(2)-pi/2)],'color',[0 0 0 0.1],'LineWidth',3);
        % cart
        xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'k','FaceAlpha',0.1)

        % ground
        xp = [x_min x_max x_max x_min ];
        yp = [-2 -2 0 0];
        patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');

        axis equal
        xlim([x_min x_max])
        ylim([-1 2])
        text(-1.5,1.5,['Time: ' num2str(t_grid(ii),'%.2f') ' s'],'interpreter','latex','fontsize',15)

        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount', inf,'DelayTime', time_step);
        else
            imwrite(imind,cm,filename,'gif', 'WriteMode', 'append','DelayTime', time_step);
        end

        if ii<length(q1_opt)
            clf;
        end
    end

    %% states
    figure
    subplot(311)
    plot(t_grid,q1_opt)
    hold on
    plot(t_grid,q2_opt)
    yline(pi,'k--');
    ylabel('$q(t)$','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    grid on
    legend({'$q_1(t)$ - cart ','$q_2(t)$ - pole'},'Interpreter','latex','Location','best')
    subplot(312)
    plot(t_grid,v1_opt)
    hold on
    plot(t_grid,v2_opt)
    yline(0,'k--');
    ylabel('$v(t)$','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    grid on
    legend({'$v_1(t)$ - cart ','$v_2(t)$ - pole'},'Interpreter','latex','Location','best')
    subplot(313)
    stairs(t_grid_u,[u_opt,nan]);
    ylabel('$u(t)$','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    grid on
end
