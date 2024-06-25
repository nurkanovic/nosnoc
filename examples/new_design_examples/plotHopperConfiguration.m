function plotHandle = plotHopperConfiguration(varargin)
q = varargin{1};
if nargin > 1
    alpha = varargin{2};
else
    alpha = 1;
end
x_head = q(1);
y_head = q(2);
y_foot = q(2) - q(4).*cos(q(3));
x_foot = q(1) - q(4).*sin(q(3));
robot_color = [0 0 0];
  % The Robot
    % head
    h = plot(x_head,y_head,'Color',[robot_color alpha],'Marker','o','MarkerSize',10,'MarkerFaceColor',(1-alpha)+robot_color,'MarkerEdgeColor',(1-alpha)+robot_color);
    hold on
    % foot
    plot(x_foot,y_foot,'Color',[robot_color alpha],'Marker','o','MarkerSize',4,'MarkerFaceColor',(1-alpha)+robot_color,'MarkerEdgeColor',(1-alpha)+robot_color);
    % link
    plotHandle = plot([x_head,x_foot],[y_head,y_foot],'Color',[robot_color alpha],'LineWidth',3);
end