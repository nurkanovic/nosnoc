%% benchmark_globals

% dream values
NS_VALUES = [1, 2, 3, 4];
NSIM_VALUES = [80, 160, 320];
NFE_VALUES = [2, 4, 8, 10, 12, 20];
IRK_SCHEME = IRKSchemes.GAUSS_LEGENDRE;

% test values
% NS_VALUES = [2];
% NSIM_VALUES = [160];
% NFE_VALUES = [2];
% IRK_SCHEME = IRKSchemes.GAUSS_LEGENDRE;

T_sim = 1;
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
e = 0.8;

x0 = [1;2;0;0];