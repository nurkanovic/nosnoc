function plot_results_ms(model, problem_options, ocp_solver, q_target, path_constraint, track_width, omega, chicane_tightness, chicane_width)
    %% read solutions
    u_opt = ocp_solver.get("u");
    u1_opt = u_opt(1,:);
    u2_opt = u_opt(2,:);

    x_opt = ocp_solver.get("x");
    x1_opt = x_opt(1,:);
    x2_opt = x_opt(2,:);
    x3_opt = x_opt(3,:);
    x4_opt = x_opt(4,:);
    x5_opt = x_opt(5,:);
    
    t_grid = ocp_solver.get_time_grid();
    u_grid = ocp_solver.get_control_grid();
    %% controls
    figure
    stairs(u_grid,[u1_opt,nan])
    hold on
    stairs(u_grid,[u2_opt,nan])
    xlabel('$t$','interpreter','latex');
    ylabel('$u(t)$','interpreter','latex');
    ylim([-2.2 2.2])
    grid on

    %% plots
    figure
    subplot(211)
    plot(x1_opt,x2_opt);
    hold on
    plot(x1_opt,x2_opt,'r.');
    xlabel('$q_x$','interpreter','latex');
    ylabel('$q_y$','interpreter','latex');
    hold on
    grid on
    plot(q_target(1),q_target(2),'rx','MarkerSize',6)
    xx = linspace(0,q_target(1),1e2);

    switch path_constraint
      case 'none'
      case 'linear'
        plot(xx ,(xx)-track_width,'k');
        plot(xx ,(xx)+track_width,'k');
      case 'nonlinear'
        plot(xx ,sin(omega*xx)-track_width,'k');
        plot(xx ,sin(omega*xx)+track_width,'k');
      case 'track'
        % figure
        arg1 = xx-pi;
        arg2 = xx-2*pi;
        sig = 1e-1;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        yy = sin(xx).*(1-step1)+(pi-xx).*step1.*(1-step2)+(-pi-sin(xx)).*step2;
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
      case 'chicane'
        yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
    end
    grid on
    axis equal

    subplot(212)
    stairs(x1_opt(1:problem_options.N_finite_elements(1):end),[u1_opt,nan]);
    hold on
    stairs(x1_opt(1:problem_options.N_finite_elements(1):end),[u2_opt,nan]);
    xlabel('$u(q_x)$','interpreter','latex');
    ylabel('$q_x$','interpreter','latex');
    legend({'$a(x)$','$s(x)$'},'interpreter','latex');
    ylim([-2.2 2.2])
    grid on


    %%  Normal Force, Veolocity,
    v_tangent = [];
    v_normal  = [];
    v_opt = [x3_opt;x4_opt];
    tangent = [cos(x5_opt);sin(x5_opt)];
    normal = [-sin(x5_opt);cos(x5_opt)];
    for ii = 1 :length(v_opt);
        v_tangent = [v_tangent,tangent(:,ii)'*v_opt(:,ii)];
        v_normal  = [v_normal,normal(:,ii)'*v_opt(:,ii)];
    end
    figure
    plot(t_grid,v_normal);
    hold on
    plot(t_grid,v_tangent);
    xlabel('$t$','interpreter','latex');
    ylabel('$v$','interpreter','latex');
    legend({'$n^\top v$','$t^\top v$'},'interpreter','latex');

end

