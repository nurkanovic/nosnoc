function model = unitri_ai_model() 
%% quadroped model equations
import casadi.*
%% integrator setings
% model
% dimensoon
n_dim = 2;
% q = [x; z; torso; thigh_1; calf_1; thigh_2; calf_2; thigh_3; calf_3; thigh_4; calf_4] 
% Dimensions
n_q = 2 + 5 + 4;            % configuration dimension
n_u = 4 + 4;                % control dimension
n_c = 4;                    % number of contact points
n_f = 2;                   % number of parameters for friction cone
% n_b = nc * nf;
% n_s = 1;

% start poulating the model;
model.n_q = n_q;
model.n_dim = n_dim;
model.n_u = n_u;

q = SX.sym('q',n_q);
v = SX.sym('v',n_q);
u = SX.sym('u',n_u);
% start poulating the model;
model.q = q;
model.v = v;
model.u = u;
model.x = [q;v];

q0 = [0.0000
    0.3236
    1.5708
   -0.6283
    0.6283
   -0.6283
    0.6283
   -0.6283
    0.6283
   -0.6283
    0.6283];

v0 = zeros(n_q,1);

x0 = [q0;v0];


model.q0 = q0;
model.v0 = v0;
model.x0 = x0;
% generalized position
model.q_x = q(1);
model.q_z = q(2);
model.q_torso = q(3);
model.q_thigh_1 = q(4);
model.q_calf_1 = q(5);
model.q_thigh_2 = q(6);
model.q_calf_2 = q(7);
model.q_thigh_3 = q(8);
model.q_calf_3 = q(9);
model.q_thigh_4 = q(10);
model.q_calf_4 = q(11);
% generalized velocity
model.v_x = v(1);
model.v_z = v(2);
model.v_torso = v(3);
model.v_thigh_1 = v(4);
model.v_calf_1 = v(5);
model.v_thigh_2 = v(6);
model.v_calf_2 = v(7);
model.v_thigh_3 = v(8);
model.v_calf_3 = v(9);
model.v_thigh_4 = v(10);
model.v_calf_4 = v(11);
% --8 control input(4 leg, each leg has two input: thigh and cart) 
model.tau_1 = u(1); % tau1: control input for thigh 1
model.tau_2 = u(2); % tau2: control input for calf 1
model.tau_3 = u(3); % tau3: control input for thigh 2
model.tau_4 = u(4);% tau4: control input for calf 2
model.tau_5 = u(5); % tau5: control input for thigh 3
model.tau_6 = u(6);%  tau 6: control input for calf 3
model.tau_7 = u(7); % tau7: control input for thigh 4
model.tau_8 = u(8); % tau8: control input for calf 4
model.ubu = 33.5*ones(n_u,1);
model.lbu = -33.5*ones(n_u,1);


model.ubx =  [1; 1.2; 1*pi; 1*pi*ones(8, 1);50 * ones(11, 1)];  
model.lbx = [0; 0.2; -1*pi; -1*pi*ones(8, 1);-50 * ones(11, 1)]; 

% World parameters
e = 0; % all inelastic impcats
mu = 0.5;      % coefficient of friction
g = 9.81*0;     % gravity
jointFriction = 0.1;

% puplate the model
model.e = e;
model.mu = mu;
model.g = g;
model.jointFriction = jointFriction;

%% model A1
% Model parameters
% Mass
m_torso = 4.713 + 4 * 0.696;
m_thigh = 1.013;
m_calf = 0.166;
mass = [m_torso;m_thigh;m_calf];
model.mass = mass;
% Interita
I_torso = 0.056579028 + 4 * 0.696 * 0.183^2.0;
I_thigh = 0.005139339;
I_calf = 0.003014022;
inertia = [I_torso;I_thigh;I_calf];
model.inertia = inertia;

% linkLength
l_torso = 0.267;
l_thigh = 0.2;
l_leg = 0.2;
linkLength = [l_torso;l_thigh;l_leg];
model.linkLength = linkLength;


% linkCenter  
d_torso = 0.5 * l_torso + 0.0127;
d_thigh = 0.5 * l_thigh - 0.00323;
d_leg = 0.5 * l_leg - 0.006435;
linkCenter = [d_torso;d_thigh;d_leg];
model.linkCenter = linkLength;

%% Inertia matrix
M = diag([0, 0, I_torso, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf]);
% torso
J_torso = jacobian_1(model, q, 'torso', 'com');
M = M + m_torso * J_torso' * J_torso;

% thigh 1
J_thigh_1 = jacobian_1(model, q, 'thigh_1', 'com');
M = M + m_thigh * J_thigh_1' * J_thigh_1;

% leg 1
J_calf_1 = jacobian_2(model, q, 'calf_1', 'com');
M = M + m_calf * J_calf_1' * J_calf_1;

% thigh 2
J_thigh_2 = jacobian_1(model, q, 'thigh_2', 'com');
M = M + m_thigh * J_thigh_2' * J_thigh_2;

% leg 2
J_calf_2 = jacobian_2(model, q, 'calf_2', 'com');
M = M + m_calf * J_calf_2' * J_calf_2;

% thigh 3
J_thigh_3 = jacobian_2(model, q, 'thigh_3', 'com');
M = M + m_thigh * J_thigh_3' * J_thigh_3;

% leg 3
J_calf_3 = jacobian_3(model, q, 'calf_3', 'com');
M = M + m_calf * J_calf_3' * J_calf_3;

% thigh 4
J_thigh_4 = jacobian_2(model, q, 'thigh_4', 'com');
M = M + m_thigh * J_thigh_4' * J_thigh_4;

% leg 4
J_calf_4 = jacobian_3(model, q, 'calf_4', 'com');
M = M + m_calf * J_calf_4' * J_calf_4;

model.M = M;
%% Forces
C = C_func(model, q, v);
B = B_func();
% H = -C + B' * tau_control + W_N' * [pN_1; pN_2; pN_3; pN_4] + W_T' * [pT_1; pT_2; pT_3; pT_4];

% disp('ignoring controls for now....')
f_v = -C+B'*u;
% f_v = -C;
%% Switching functions
f_c = gap_func(model, q)';
% Jacobians
J_normal = W_N_func(model, q);
J_tangent = W_T_func(model, q);

%% Populate model
model.f_v = f_v;
model.f_c = f_c;
model.J_normal  = J_normal';
model.J_tangent = J_tangent';

% save('unitree_ai','model');

%% kinemaics and Jacobian
function kine = kinematics_1(model, q, body, mode)
x = q(1);
z = q(2);
if strcmp(body, 'torso')
    l = model.linkLength(1);
    d = model.linkCenter(1);
elseif strcmp(body, 'thigh_1') || strcmp(body, 'thigh_2')
    l = model.linkLength(2);
    d = model.linkCenter(2);
end

switch body
    case 'torso'
        q_arg = q(3);        
    case 'thigh_1'
        q_arg = q(4);        
    case 'thigh_2'
        q_arg = q(6);               
end

switch mode
    case 'ee'
        kine = [x + l * sin(q_arg); z - l * cos(q_arg)];
    case 'com'
        kine = [x + d * sin(q_arg); z - d * cos(q_arg)];
end

end

function jac = jacobian_1(model, q, body, mode)
if strcmp(body, 'torso')
    l = model.linkLength(1);
    d = model.linkCenter(1);
elseif strcmp(body, 'thigh_1') || strcmp(body, 'thigh_2')
    l = model.linkLength(2);
    d = model.linkCenter(2);
end

switch body
    case 'torso'
        q_arg = q(3);        
        switch mode
            case 'ee'
                jac = [1, 0, l * cos(q_arg), 0, 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, l * sin(q_arg), 0, 0, 0, 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, d * cos(q_arg), 0, 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, d * sin(q_arg), 0, 0, 0, 0, 0, 0, 0, 0];                
        end
    case 'thigh_1'
        q_arg = q(4); 
        switch mode
            case 'ee'
                jac = [1, 0, 0, l * cos(q_arg), 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, 0, l * sin(q_arg), 0, 0, 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, 0, d * cos(q_arg), 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, 0, d * sin(q_arg), 0, 0, 0, 0, 0, 0, 0];                
        end
    case 'thigh_2'
        q_arg = q(6);   
        switch mode
            case 'ee'
                jac = [1, 0, 0, 0, 0, l * cos(q_arg), 0, 0, 0, 0, 0;...
                       0, 1, 0, 0, 0, l * sin(q_arg), 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, 0, 0, 0, d * cos(q_arg), 0, 0, 0, 0, 0;...
                       0, 1, 0, 0, 0, d * sin(q_arg), 0, 0, 0, 0, 0];                
        end
end

end

function kine = kinematics_2(model, q, body, mode)
if strcmp(body, 'calf_1') || strcmp(body, 'calf_2')
    l = model.linkLength(3);
    d = model.linkCenter(3);
elseif strcmp(body, 'thigh_3') || strcmp(body, 'thigh_4')
    l = model.linkLength(2);
    d = model.linkCenter(2);
end

switch body
    case 'calf_1'
        p = kinematics_1(model, q, 'thigh_1', 'ee');
        q_arg = q(5);
    case 'calf_2'
        p = kinematics_1(model, q, 'thigh_2', 'ee');
        q_arg = q(7);      
    case 'thigh_3'
        p = kinematics_1(model, q, 'torso', 'ee');
        q_arg = q(8);
    case 'thigh_4'
        p = kinematics_1(model, q, 'torso', 'ee');
        q_arg = q(10);
end

switch mode
    case 'ee'
        kine = p + [l * sin(q_arg); -l*cos(q_arg)];
    case 'com'
        kine = p + [d * sin(q_arg); -d*cos(q_arg)];
end

end

function jac = jacobian_2(model, q, body, mode)
if strcmp(body, 'calf_1') || strcmp(body, 'calf_2')
    l = model.linkLength(3);
    d = model.linkCenter(3);
elseif strcmp(body, 'thigh_3') || strcmp(body, 'thigh_4')
    l = model.linkLength(2);
    d = model.linkCenter(2);
end

switch body
    case 'calf_1'
        jac = jacobian_1(model, q, 'thigh_1', 'ee');
        q_arg = q(5);
        switch mode
            case 'ee'
                jac(1, 5) = jac(1, 5) + l * cos(q_arg);
                jac(2, 5) = jac(2, 5) + l * sin(q_arg);
            case 'com'
                jac(1, 5) = jac(1, 5) + d * cos(q_arg);
                jac(2, 5) = jac(2, 5) + d * sin(q_arg);                
        end
    case 'calf_2'
        jac = jacobian_1(model, q, 'thigh_2', 'ee');
        q_arg = q(7);
        switch mode
            case 'ee'
                jac(1, 7) = jac(1, 7) + l * cos(q_arg);
                jac(2, 7) = jac(2, 7) + l * sin(q_arg);
            case 'com'
                jac(1, 7) = jac(1, 7) + d * cos(q_arg);
                jac(2, 7) = jac(2, 7) + d * sin(q_arg);                
        end       
    case 'thigh_3'
        jac = jacobian_1(model, q, 'torso', 'ee');
        q_arg = q(8);
        switch mode
            case 'ee'
                jac(1, 8) = jac(1, 8) + l * cos(q_arg);
                jac(2, 8) = jac(2, 8) + l * sin(q_arg);
            case 'com'
                jac(1, 8) = jac(1, 8) + d * cos(q_arg);
                jac(2, 8) = jac(2, 8) + d * sin(q_arg);                
        end
    case 'thigh_4'
        jac = jacobian_1(model, q, 'torso', 'ee');
        q_arg = q(10);
        switch mode
            case 'ee'
                jac(1, 10) = jac(1, 10) + l * cos(q_arg);
                jac(2, 10) = jac(2, 10) + l * sin(q_arg);
            case 'com'
                jac(1, 10) = jac(1, 10) + l * cos(q_arg);
                jac(2, 10) = jac(2, 10) + l * sin(q_arg);
        end
end

end

function kine = kinematics_3(model, q, body, mode)
l = model.linkLength(3);
d = model.linkCenter(3);
switch body
    case 'calf_3'
        p = kinematics_2(model, q, 'thigh_3', 'ee');
        q_arg = q(9);
    case 'calf_4'
        p = kinematics_2(model, q, 'thigh_4', 'ee');
        q_arg = q(11);
end

switch mode
    case 'ee'
        kine = p + [l * sin(q_arg); -l * cos(q_arg)];
    case 'com'
        kine = p + [d * sin(q_arg); -d * cos(q_arg)];
end

end

function jac = jacobian_3(model, q, body, mode)
l = model.linkLength(3);
d = model.linkCenter(3);
switch body
    case 'calf_3'
        jac = jacobian_2(model, q, 'thigh_3', 'ee');
        q_arg = q(9);
        switch mode
            case 'ee'
                jac(1, 9) = jac(1, 9) + l * cos(q_arg);
                jac(2, 9) = jac(2, 9) + l * sin(q_arg);
            case 'com'
                jac(1, 9) = jac(1, 9) + d * cos(q_arg);
                jac(2, 9) = jac(2, 9) + d * sin(q_arg);
        end
    case 'calf_4'
        jac = jacobian_2(model, q, 'thigh_4', 'ee');
        q_arg = q(11);
        switch mode
            case 'ee'
                jac(1, 11) = jac(1, 11) + l * cos(q_arg);
                jac(2, 11) = jac(2, 11) + l * sin(q_arg);
            case 'com'
                jac(1, 11) = jac(1, 11) + d * cos(q_arg);
                jac(2, 11) = jac(2, 11) + d * sin(q_arg);
        end
end

end

%% Lagrangian
function L = lagrangian(model, q, v)
m_torso = model.mass(1);
m_thigh = model.mass(2);
m_calf = model.mass(3);
I_torso = model.inertia(1);
I_thigh = model.inertia(2);
I_calf = model.inertia(3);
g = model.g;

L = 0;

% torso
p_torso = kinematics_1(model, q, 'torso', 'com');
J_torso = jacobian_1(model, q, 'torso', 'com');
v_torso = J_torso * v;

L = L + 0.5 * m_torso * v_torso' * v_torso;
L = L + 0.5 * I_torso * v(3)^2;
L = L - m_torso * g * p_torso(2);

% thigh 1
p_thigh_1 = kinematics_1(model, q, 'thigh_1', 'com');
J_thigh_1 = jacobian_1(model, q, 'thigh_1', 'com');
v_thigh_1 = J_thigh_1 * v;

L = L + 0.5 * m_thigh * v_thigh_1' * v_thigh_1;
L = L + 0.5 * I_thigh * v(4)^2;
L = L - m_thigh * g * p_thigh_1(2);

% leg 1
p_calf_1 = kinematics_2(model, q, 'calf_1', 'com');
J_calf_1 = jacobian_2(model, q, 'calf_1', 'com');
v_calf_1 = J_calf_1 * v;

L = L + 0.5 * m_calf * v_calf_1' * v_calf_1;
L = L + 0.5 * I_calf * v(5)^2;
L = L - m_calf * g * p_calf_1(2);

% thigh 2
p_thigh_2 = kinematics_1(model, q, 'thigh_2', 'com');
J_thigh_2 = jacobian_1(model, q, 'thigh_2', 'com');
v_thigh_2 = J_thigh_2 * v;

L = L + 0.5 * m_thigh * v_thigh_2' * v_thigh_2;
L = L + 0.5 * I_thigh * v(6)^2;
L = L - m_thigh * g * p_thigh_2(2);

% leg 2
p_calf_2 = kinematics_2(model, q, 'calf_2', 'com');
J_calf_2 = jacobian_2(model, q, 'calf_2', 'com');
v_calf_2 = J_calf_2 * v;

L = L + 0.5 * m_calf * v_calf_2' * v_calf_2;
L = L + 0.5 * I_calf * v(7)^2;
L = L - m_calf * g * p_calf_2(2);

% thigh 3
p_thigh_3 = kinematics_2(model, q, 'thigh_3', 'com');
J_thigh_3 = jacobian_2(model, q, 'thigh_3', 'com');
v_thigh_3 = J_thigh_3 * v;

L = L + 0.5 * m_thigh * v_thigh_3' * v_thigh_3;
L = L + 0.5 * I_thigh * v(8)^2;
L = L - m_thigh * g * p_thigh_3(2);

% leg 3
p_calf_3 = kinematics_3(model, q, 'calf_3', 'com');
J_calf_3 = jacobian_3(model, q, 'calf_3', 'com');
v_calf_3 = J_calf_3 * v;

L = L + 0.5 * m_calf * v_calf_3' * v_calf_3;
L = L + 0.5 * I_calf * v(9)^2;
L = L - m_calf * g * p_calf_3(2);

% thigh 4
p_thigh_4 = kinematics_2(model, q, 'thigh_4', 'com');
J_thigh_4 = jacobian_2(model, q, 'thigh_4', 'com');
v_thigh_4 = J_thigh_4 * v;

L = L + 0.5 * m_thigh * v_thigh_4' * v_thigh_4;
L = L + 0.5 * I_thigh * v(10)^2;
L = L - m_thigh * g * p_thigh_4(2);

% leg 4
p_calf_4 = kinematics_3(model, q, 'calf_4', 'com');
J_calf_4 = jacobian_3(model, q, 'calf_4', 'com');
v_calf_4 = J_calf_4 * v;

L = L + 0.5 * m_calf * v_calf_4' * v_calf_4;
L = L + 0.5 * I_calf * v(11)^2;
L = L - m_calf * g * p_calf_4(2);

end

function dLv = dLv_func(model, q, v)

L = lagrangian(model, q, v);
dLv = jacobian(L, q);

end

function dLdv = dLdv_func(model, q, v)

L = lagrangian(model, q, v);
dLdv = jacobian(L, v);

end

%% dynamics matrix
function M = M_func(model, q)
m_torso = model.mass(1);
m_thigh = model.mass(2);
m_calf = model.mass(3);
I_torso = model.inertia(1);
I_thigh = model.inertia(2);
I_calf = model.inertia(3);

M = diag([0, 0, I_torso, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf]);
% torso
J_torso = jacobian_1(model, q, 'torso', 'com');
M = M + m_torso * J_torso' * J_torso;

% thigh 1
J_thigh_1 = jacobian_1(model, q, 'thigh_1', 'com');
M = M + m_thigh * J_thigh_1' * J_thigh_1;

% leg 1
J_calf_1 = jacobian_2(model, q, 'calf_1', 'com');
M = M + m_calf * J_calf_1' * J_calf_1;

% thigh 2
J_thigh_2 = jacobian_1(model, q, 'thigh_2', 'com');
M = M + m_thigh * J_thigh_2' * J_thigh_2;

% leg 2
J_calf_2 = jacobian_2(model, q, 'calf_2', 'com');
M = M + m_calf * J_calf_2' * J_calf_2;

% thigh 3
J_thigh_3 = jacobian_2(model, q, 'thigh_3', 'com');
M = M + m_thigh * J_thigh_3' * J_thigh_3;

% leg 3
J_calf_3 = jacobian_3(model, q, 'calf_3', 'com');
M = M + m_calf * J_calf_3' * J_calf_3;

% thigh 4
J_thigh_4 = jacobian_2(model, q, 'thigh_4', 'com');
M = M + m_thigh * J_thigh_4' * J_thigh_4;

% leg 4
J_calf_4 = jacobian_3(model, q, 'calf_4', 'com');
M = M + m_calf * J_calf_4' * J_calf_4;

end

function C = C_func(model, q, v)

dLv = dLv_func(model, q, v);
dLdv = dLdv_func(model, q, v);

dLdv_v = jacobian(dLdv, q);

C =  dLdv_v*v - dLv';
end

function gap = gap_func(model, q)
p_calf_1 = kinematics_2(model, q, 'calf_1', 'ee');
p_calf_2 = kinematics_2(model, q, 'calf_2', 'ee');
p_calf_3 = kinematics_3(model, q, 'calf_3', 'ee');
p_calf_4 = kinematics_3(model, q, 'calf_4', 'ee');

gap = [p_calf_1(2), p_calf_2(2), p_calf_3(2), p_calf_4(2)];
end

function B = B_func()
B = [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0;...
     0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0;...
     0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0;...
     0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0;...
     0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0;...
     0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0;...
     0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0;...
     0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1];
end

function W_N = W_N_func(model, q)
J_calf_1 = jacobian_2(model, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(model, q, 'calf_2', 'ee');
J_calf_3 = jacobian_3(model, q, 'calf_3', 'ee');
J_calf_4 = jacobian_3(model, q, 'calf_4', 'ee');

W_N = [J_calf_1(2, :);...
       J_calf_2(2, :);...
       J_calf_3(2, :);...
       J_calf_4(2, :)];
end

function W_T = W_T_func(model, q)
J_calf_1 = jacobian_2(model, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(model, q, 'calf_2', 'ee');
J_calf_3 = jacobian_3(model, q, 'calf_3', 'ee');
J_calf_4 = jacobian_3(model, q, 'calf_4', 'ee');

W_T = [J_calf_1(1, :);...
       J_calf_2(1, :);...
       J_calf_3(1, :);...
       J_calf_4(1, :)];
end

end