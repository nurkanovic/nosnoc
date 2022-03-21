function [model,settings] = model_reformulation_fesd(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
unfold_struct(settings,'caller')


%% Step-size
model.h = T/N_stages;
%% Check are x and u provided
if exist('x')
    n_x = length(x);
    % check  lbx
    if exist('lbx')
        if length(lbx) ~= n_x
            error('The vector lbx, for the lower bounds of x has the wrong size.')
        end
    else
        lbx = -inf*ones(n_x,1);
    end
    % check ubx
    if exist('ubx')
        if length(ubx) ~= n_x
            error('The vector ubx, for the upper bounds of x has the wrong size.')
        end
    else
        ubx = inf*ones(n_x,1);
    end

else
    error('Please provide the state vector x, a CasADi symbolic variable.')
end

if exist('u')
    n_u = length(u);
    % check  lbu
    if exist('lbu')
        if length(lbu) ~= n_u
            error('The vector lbu, for the lower bounds of u has the wrong size.')
        end
    else
        lbu = -inf*ones(n_u,1);
    end
    % check ubu
    if exist('ubu')
        if length(ubu) ~= n_u
            error('The vector ubu, for the upper bounds of u has the wrong size.')
        end
    else
        ubu = inf*ones(n_u,1);
    end

    % check u0
    if exist('u0')
        if length(u0) ~= n_u
            error('The vector u0, for the initial guess of u has the wrong size.')
        end
    else
        u0 = 0*ones(n_u,1);
    end
    %
    %     if simulation_problem
    %             lbu = 0*ones(n_u,1);
    %             ubu = 0*ones(n_u,1);
    %     end
else
    warning('No control vector u is provided.')
    lbu = [];
    ubu = [];
end
%% Stage and terminal costs
if ~exist('f_q')
    fprintf('Info: No stage cost is provided. \n')
    f_q = 0;
end

if exist('f_q_T')
    terminal_cost = 1;
else
    fprintf('Info: No terminal cost is provided. \n')
    f_q_T = 0;
end

%% Inequality constraints
if exist('g_ineq')
    general_nonlinear_constraint  = 1;
    n_g_ineq = length(g_ineq);
    if exist('g_ineq_lb')
        if length(g_ineq_lb)~=n_g_ineq;
            error('The user provided vector g_ineq_lb has the wrong size.')
        end
    else
        g_ineq_lb = -inf*ones(n_g_ineq,1);
    end

    if exist('g_ineq_ub')
        if length(g_ineq_ub)~=n_g_ineq;
            error('The user provided vector g_ineq_ub has the wrong size.')
        end
    else
        g_ineq_ub =  0*ones(n_g_ineq,1);
    end
    g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});
else
    n_g_ineq = 0;
    general_nonlinear_constraint  = 0;
    fprintf('Info: No path constraints are provided. \n')
end
%% Terminal constraints
if exist('g_terminal')
    terminal_constraint = 1;
    n_g_terminal = length(g_terminal);
    if exist('g_terminal_lb')
        if length(g_terminal_lb)~=n_g_terminal;
            error('The user provided vector g_terminal_lb has the wrong size.')
        end
    else
        g_terminal_lb = 0*ones(n_g_terminal,1);
    end

    if exist('g_terminal_ub')
        if length(g_terminal_ub)~=n_g_terminal;
            error('The user provided vector g_terminal_ub has the wrong size.')
        end
    else
        g_terminal_ub =  0*ones(n_g_terminal,1);
    end
    g_terminal_fun  = Function('g_terminal_fun',{x},{g_terminal});
else
    terminal_constraint = 0;
    n_g_terminal = 0;
    fprintf('Info: No terminal constraints are provided. \n')
end
%% Default settings for submodes , m simplex (Carteisan product of Filippov systems in FESD Paper)
if ~exist('n_simplex');
    n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
end



%% Stewart's representation of the sets R_i and discirimant functions g_i (here called h_i)

if ~exist('F')
    error('Matrix F with PSS modes not provided.');
else
    m_1 = size(F,2);
    m_vec = [m_1];
end

if ~exist('S')
    if exist('g_discriminant')
        g_discriminant_1 = g_discriminant;
    else
        error('Sign matrix S for regions is not provided.');
    end
else
    if ~exist('c')
        error('Region constraint function not provided.');
    end

    if size(S,2) ~= length(c)
        error('The matrix S and vector c do not have compatible dimension.');
    end
    % discrimnant functions
    g_discriminant_1 = -S*c;

end

h_1 = g_discriminant_1; % TODO: fix notation
g_discriminants = [h_1];
f_1 = F;

%% Declare model variables and equations
m_ind_vec = [cumsum(m_vec)-m_1+1]; % index ranges of the corresponding thetas and lambdas
m = sum(m_vec);
%% Parameters
sigma = MX.sym('sigma'); % homotopy parameter
n_param = 1;  % number of parameters,  we model it as control variables and merge them with simple equality constraints
p = [sigma];
n_p = 1;

%% Algebraic variables defintion

n_theta = sum(m_vec); % number of modes
n_lambda = n_theta;
n_z = n_theta+n_lambda+n_simplex; % n_theta + n_lambda + n_mu

theta = [];
mu = [];
lambda = [];

E = []; % diagonal matrix with one vectors

for i = 1:n_simplex
    i_str = num2str(i);
    % define theta (convex multiplers of Filippov embedding)
    eval(['theta_' i_str '= MX.sym(''theta_' i_str ''',m_' i_str ');'])
    eval(['theta = [theta; theta_' i_str '];'])
    % define mu_i (Lagrange multipler of e'theta =1;)
    eval(['mu_' i_str '= MX.sym(''mu_' i_str ''',1);'])
    eval(['mu = [mu; mu_' i_str '];'])
    % define lambda_i (Lagrange multipler of theta >= 0;)
    eval(['lambda_' i_str '= MX.sym(''lambda_' i_str ''',m_' i_str ');'])
    eval(['lambda = [lambda; lambda_' i_str '];'])
    % adefine ppropiate vector of ones
    eval(['e_' i_str '=ones(m_' i_str ' ,1);'])

    if n_simplex > 1
        eval(['E_row = zeros(m_' i_str ',n_simplex);'])
        eval(['E_row(:,' i_str ') = e_' i_str ';']);
        E = [E;E_row];
    end

end
if n_simplex == 1
    E = e_1;
end

% define algerbraic variables which arise from Stewart's reformulation of a PSS into a DCS
z = [theta;lambda;mu];

%% Reformulate the Filippov ODE into a DCS
f_x = zeros(n_x,1);
% rhs of ODE;
for i = 1:n_simplex
    i_str = num2str(i);
    try
        eval(['f_x = f_x + f_' i_str '*theta_'  i_str  ';']);
    catch
        eval(['f_x = f_x+ F*theta_'  i_str  ';']);
    end
end

% basic algebraic equations and complementarty condtions of the DCS
% (Note that the cross complementarities are later defined when the discrete
% time variables for every IRK stage in the create_nlp_fesd function are defined.)

% h_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_simplex
% lambda_i'*theta_i = 0; for all i = 1,..., n_simplex
% lambda_i >= 0;    for all i = 1,..., n_simplex
% theta_i >= 0;     for all i = 1,..., n_simplex

f_z = []; % collects standard algebraic equations 0 = g_i(x) - \lambda_i - e \mu_i
f_z_convex = []; % equation for the convex multiplers 1 = e' \theta
f_z_cc = []; % the orthogonality conditions diag(\theta) \lambda = 0.
for i = 1:n_simplex
    i_str = num2str(i);
    % Gradient of Lagrange Function of indicator LP
    eval(['f_z = [f_z; h_' i_str '-lambda_'  i_str '+mu_'  i_str '*e_' i_str '];']);
    % Sum of all thethas equal to 1
    eval(['f_z_convex = [f_z_convex; ; e_' i_str '''*theta_'  i_str '-1];']);
    % Complementarty conditions arising from theta_i >= 0
    %     if pointwise_or_integral
    eval(['f_z_cc = [f_z_cc; lambda_' i_str '.*theta_'  i_str '];']);
    %     else
    %         eval(['f_z_cc = [f_z_cc; lambda_' i_str '''*theta_'  i_str '];']);
    %     end
end

% such orderdng to have better sparsity. Collect the first two equations,
% as they do not need to be relaxed or smoothed. The evaluations of f_z_cc have latter,
% depending on the MPCC algrotim a special treatment in create_nlp_fesd().
f_z = [f_z;f_z_convex];

%% MPCC Specific Considerations
if mpcc_mode == 4
    % f_z  already defined
    n_bilinear_cc = 0;
else
    %     if use_fesd
    %         %         f_z = [f_z];
    %         % all complementarites are treated in the create_nlp function since
    %         % there is a coupling of variables across different time steps.
    %     else
    %         % TODO: all at once, three, or full point_wise
    %         % complementarity constraints (if mfe off)
    %         f_z_cc = f_z_cc;
    %         %         f_z = [f_z;];
    %     end
    %     if pointwise_or_integral
    if use_fesd
        % TODO remove this,
        switch cross_complementarity_mode
            case 1
                n_bilinear_cc = n_theta;
            case 2
                n_bilinear_cc = n_theta;
            case 3
                n_bilinear_cc = (d+1)*n_theta;
            case 4
                n_bilinear_cc = (d+1)*n_theta;
            otherwise
                warning('pick fesd_complementartiy_mode between 1 and 4, setting to default value = 1');
                fesd_complementartiy_mode = 1;
                n_bilinear_cc = n_theta;
        end
    else
        n_bilinear_cc = n_theta;
    end
    %     else
    %         n_bilinear_cc = n_simplex;
    %     end
end
% sum over all complementariteies
J_cc = sum(f_z_cc);  % (used in l1 penalties and for evaluation of resiudal)
% Point-wise
n_algebraic_constraints = n_theta+n_simplex;  % dim(g  - lambda - mu *e ) + dim( E theta) ;
n_algebraic_constraints_total =n_algebraic_constraints;
% n_algebraic_constraints_total =n_algebraic_constraints+n_bilinear_cc; %

%% CasADi functions
% model equations
try
    h_fun = Function('h_fun',{x},{g_discriminants});
    c_fun = Function('c_fun',{x},{c});
catch
    h_fun = Function('h_fun',{x,u},{g_discriminants});
    c_fun = Function('c_fun',{x,u},{c});
end
%% TODO: fix with fun and without fun --> make consistsnet
f_x = Function('f_x',{x,z,u},{f_x,f_q});
f_z = Function('f_z',{x,z,u},{f_z});
f_z_cc = Function('f_z_cc',{x,z,u},{f_z_cc});
f_q_T = Function('f_q_T',{x},{f_q_T});
% l1 norm of all complementarity pairs;
f_J_cc = Function('f_J_cc',{x,z,u},{J_cc});
%% Function for continious time output of algebraic variables.
[B,C,D,tau_root] = collocation_times_fesd(d,collocation_scheme);
tau_col = tau_root(2:end);
tau = MX.sym('tau');  % Time argument
h_scale = MX.sym('h_scale');  % Rescale the intevral from [0 1] to [0 h_scale]; h_scale will usually be the step size
y = MX.sym('y',d);
y_vec = MX.sym('y_vec',n_theta,d);
lagrange_poly = 0;
m_lagrange = [];
for ii = 1:d
    term_ii = y(ii);
    m_temp = 1;
    tau_scaled= h_scale*tau;
    tau_col_scaled = h_scale*tau_col;
    for jj = 1:d
        if jj~=ii
            term_ii = term_ii.*((tau_scaled-tau_col_scaled(jj))./(tau_col_scaled(ii)-tau_col_scaled(jj)));
            m_temp  = m_temp*((1-tau_col(jj))./(tau_col(ii)-tau_col(jj)));
        end
    end
    m_lagrange = [m_lagrange,m_temp];
    lagrange_poly = lagrange_poly + term_ii;
end

lagrange_poly_fun = Function('lagrange_poly',{y,tau,h_scale},{lagrange_poly}); % note that this is a scalar function
% lagrange_poly_fun([3 2 1],1,1)

%% Collect Outputs

%
model.sigma = sigma;
model.p = p;

model.lbx = lbx;
model.ubx = ubx;

model.lbu = lbu;
model.ubu = ubu;
model.u0 = u0;

if general_nonlinear_constraint
    model.g_ineq_lb = g_ineq_lb;
    model.g_ineq_ub = g_ineq_ub;
    model.g_ineq_fun = g_ineq_fun;
    model.general_nonlinear_constraint = general_nonlinear_constraint;
end

if terminal_constraint
    model.g_terminal_lb = g_terminal_lb;
    model.g_terminal_ub = g_terminal_ub;
    model.terminal_constraint = terminal_constraint;
    model.g_terminal_fun = g_terminal_fun;
end


model.f_x = f_x;
model.f_z = f_z;
model.f_q_T = f_q_T;
model.f_z_cc = f_z_cc;
model.f_J_cc = f_J_cc;
model.h_fun = h_fun;
model.c_fun = c_fun;

model.lagrange_poly_fun = lagrange_poly_fun;
model.m_lagrange = m_lagrange;

% Some Dimensions;
model.n_x = n_x;
model.n_z = n_z;
model.n_u = n_u;
model.n_p = n_p;
model.n_simplex = n_simplex;



model.z = z;
model.E = E;

model.m_vec = m_vec;
model.m_ind_vec = m_ind_vec;
model.n_theta = n_theta;
model.n_lambda = n_lambda;
model.n_algebraic_constraints = n_algebraic_constraints;
model.n_algebraic_constraints_total = n_algebraic_constraints_total;
model.n_bilinear_cc = n_bilinear_cc;
end
