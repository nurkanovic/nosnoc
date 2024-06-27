function plotHandle = plotHopperScene1(varargin)
x_min = varargin{1};
x_max = varargin{2};
y_min = varargin{3};
y_max = varargin{4};

if nargin > 4
    alpha = varargin{5};
else
    alpha = 0.1;
end
%% Plot
   % ground
    xp = [x_min x_max x_max x_min ];
    yp = [-1 -1 0 0];
    plotHandle = patch(xp,yp,0*ones(1,3),'FaceAlpha',alpha,'EdgeColor','none');
    axis equal
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('$x$ [m]','Interpreter','latex')
    ylabel('$y$ [m]','Interpreter','latex')

end