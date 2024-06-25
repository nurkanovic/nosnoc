function plotHandle = plotHopperSwitchingFun(varargin)
t_opt = varargin{1};
c_eval = varargin{2};
if nargin > 2
    fig_num = varargin{3};
else
    fig_num = 3;
end

figure(fig_num)
movegui('northeast');
subplot(311)
plot(t_opt,c_eval(1,:));
xlabel('$t$','interpreter','latex');
ylabel('$f_c(q)$ - gap function ','interpreter','latex');
grid on

subplot(312)
plot(t_opt,c_eval(2,:));
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_{\mathrm{n}}$ - normal velocity','interpreter','latex');

subplot(313)
plot(t_opt,c_eval(3,:));
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_{\mathrm{t}}$ -tangential velocitiy','interpreter','latex');





end