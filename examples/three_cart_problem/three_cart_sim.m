%% Three crat problem - simulation

clear all;
close all;
clc;
import casadi.*

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;
k1 = 1*0;
k2 = 1*0;
k3 = 1*0;
c_damping = 2;
u_max = 20;
u_min = -20;
ubx = [10; 15; 15; 5; 5; 5]; 
lbx = [-15; -15; -10; -5; -5; -5];            
          
x0 = [-3; 0; 3; 0; 0; 0];

%%
problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 1;
solver_options.mpcc_mode = 'elastic_ineq';
solver_options.opts_casadi_nlp.ipopt.max_iter = 5e2;
solver_options.print_level = 2;
solver_options.N_homotopy = 12;
problem_options.cross_comp_mode = 7;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 1;

problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 0;

% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%%
% Symbolic variables and bounds
q = SX.sym('q',3); v = SX.sym('v',3);
model = NosnocModel();
model.x = [q;v]; 
model.e = 0;
model.mu_f = 0;
model.a_n = 100;
model.x0 = x0; 
% fixed control
u = [-1;10;0];
model.f_v = [1/m1*(u(1)-c_damping*v(1)-k1*q(1));...
           1/m2*(u(2)-c_damping*v(2)+k2*(q(1)-q(2)));...
           1/m3*(u(3)-c_damping*v(3)+k3*q(2))];

model.f_c = [q(2) - q(1) - 0.5*cart_width2 - 0.5*cart_width1;...
           q(3) - q(2) - 0.5*cart_width3 - 0.5*cart_width2];

model.dims.n_dim_contact = 2;
%% Simulation settings
N_FE = 2;
T_sim = 3;
N_sim = 60;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();
%% read and plot results
unfold_struct(results,'base');
p1 = results.x(1,:);
p2 = results.x(2,:);
p3 = results.x(3,:);
v1 = results.x(4,:);
v2 = results.x(5,:);
v3 = results.x(6,:);
t_opt = results.x(7,:);
%%
figure
subplot(121)
plot(t_opt,p1,'LineWidth',1.5);
hold on
plot(t_opt,p2,'LineWidth',1.5);
plot(t_opt,p3,'LineWidth',1.5);
% axis equal
grid on
legend({'$p_1(t)$','$p_2(t)$','$p_3(t)$'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$p$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,v1,'LineWidth',1.5);
hold on
plot(t_opt,v2,'LineWidth',1.5);
plot(t_opt,v3,'LineWidth',1.5);
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

%% animation
figure('Renderer', 'painters', 'Position', [100 100 1400 600])

carts_appart = 2;
x_min =min([p1,p2,p3])-2.5;
x_max = max([p1,p2,p3])+2.5;
cart_height = 2;

carts_appart = 1.5*1;
for ii = 1:length(p1)
    % cart 1
    xp = [p1(ii)-cart_width1/2 p1(ii)+cart_height/2 p1(ii)+cart_height/2 p1(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.4)
    hold on
    % cart 2
    xp = [p2(ii)-cart_width2/2 p2(ii)+cart_height/2 p2(ii)+cart_height/2 p2(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.4)

    if k2 >0
    spring1 = linspace(min(p1(ii),p2(ii)),max(p1(ii),p2(ii)),5);
    plot(spring1,spring1*0+cart_height/2,'ko-','LineWidth',1);
    end

    % cart 3
    xp = [p3(ii)-cart_width3/2 p3(ii)+cart_height/2 p3(ii)+cart_height/2 p3(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.4)

    % ground     
    xp = [x_min x_max x_max x_min ];
    yp = [-1 -1 0 0];
    patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');

    
    axis equal
    xlim([x_min x_max])
    ylim([-0.75 2.5])
    pause(problem_options.h_k);
    clf
end

