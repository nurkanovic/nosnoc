function cls_callback(mpcc, nlp, opts, ipopt_solver, results)
    persistent fig_GH
    persistent GH_plot_handle
    persistent bound_plot_handle
    persistent G
    persistent H
    if isempty(fig_GH)
        fig_GH = figure();
        G = [0];
        H = [0];
        hold on
        GH_plot_handle = plot(G, H, '-o');
        GH_plot_handle.XDataSource = 'G';
        GH_plot_handle.YDataSource = 'H';
        bound_plot_handle = plot(G, H, '-r');
        bound_plot_handle.XDataSource = 'bx';
        bound_plot_handle.YDataSource = 'by';
        plot([0 10], [0 0], '-k');
        plot([0 0], [-10 10], '-k');
        hold off
        G = [];
        H = [];
        xlabel('G(w,p)');
        ylabel('H(w,p)');
    end
    % populate mpcc w
    mpcc.w.res = nlp.w.mpcc_w().res;
    % Handle GH plot
    G = [G, mpcc.w.Lambda_normal(1,2).res];
    H = [H, -mpcc.w.N_vn(1,2).res];

    x_prev = mpcc.w.x(1,1,2).res;
    x = mpcc.w.x(1,2,0).res;
    h = mpcc.w.h(1,2).res;
    bx = [-1,10];
    by = [x_prev(2)-x(1)/(0.1*h), x_prev(2)-x(1)/(0.1*h)];
    
    
    ax = gca(fig_GH);
    xlim(ax, [-1, 10]);
    ylim(ax, [-10, 10]);
    refreshdata(GH_plot_handle,'caller');
    refreshdata(bound_plot_handle,'caller');
end
