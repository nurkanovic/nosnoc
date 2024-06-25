function plotHandle = plotHopperStatesControls(varargin)
x_opt = varargin{1};
u_opt = varargin{2};
x_end = varargin{3};
if nargin > 3
    fig_num  = varargin{4};
else
    fig_num  = 7;
end


q_opt = x_opt(1:4,:);
v_opt = x_opt(5:8,:);
t_opt = x_opt(9,:);

q1 = q_opt(1,:);
q2 = q_opt(2,:);
q3 = q_opt(3,:);
q4 = q_opt(4,:);

v1 = v_opt(1,:);
v2 = v_opt(2,:);
v3 = v_opt(3,:);
v4 = v_opt(4,:);
%%
figure(fig_num)
movegui('northwest');
clf;
subplot(231)
plot(t_opt,q1)
hold on
plot(t_opt,q2)
% yline(x_end(1),'b--')
% yline(x_end(2),'r--')
xlabel('$t$','interpreter','latex');
ylabel('$q$','interpreter','latex');
legend({'$x_{\mathrm{b}}$','$y_{\mathrm{b}}$'},'interpreter','latex','location','best');
grid on
subplot(232)
plot(t_opt,q3)
% yline(x_end(3),'b--')

grid on
xlabel('$t$','interpreter','latex');
ylabel('$\theta_{\mathrm{b}}$','interpreter','latex');
grid on
subplot(233)
plot(t_opt,q4)
% yline(x_end(4),'b--')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$L_{\mathrm{link}}$','interpreter','latex');
grid on
subplot(234)
plot(t_opt,v1)
hold on
plot(t_opt,v2)
% yline(x_end(5),'b--')
% yline(x_end(6),'r--')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
legend({'$v_{x,\mathrm{b}}$','$v_{y,\mathrm{b}}$'},'interpreter','latex','location','best');
subplot(235)
plot(t_opt,v3)
% yline(x_end(7),'b--')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\omega_{\mathrm{b}}$','interpreter','latex');
% yline(x_end(8))
subplot(236)
plot(t_opt,v4)
grid on
xlabel('$t$','interpreter','latex');
ylabel('$vL_{\mathrm{link}}$','interpreter','latex');
grid on

%%  controls
N_FE = 3; % TODO: make more cleanear with model dimensions
figure(fig_num+1);
clf;
movegui('southwest');
stairs(t_opt(1:N_FE:end),[u_opt(1,:),nan])
hold on
stairs(t_opt(1:N_FE:end),[u_opt(2,:),nan])
grid on
legend({'Orientation','Horizontal'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$\tau$','interpreter','latex');
grid on


end
