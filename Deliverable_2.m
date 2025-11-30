clc; clear all; close all; 
syms theta theta_dot theta_ddot 
syms x x_dot x_ddot
%syms phi phi_dot phi_ddot
syms L0 L1 m0 m1 
syms d_phi d_theta x0 theta0 q0 FA FW theta_ref Lx
%syms g L2 m2 m3 dW dx kx k_theta

%% Parameter values
g = 9.81;
L2 = 10;
m2 = 2000;
m3 = 500;
dW = 2000;
dx = 200;
kx = 7000;
k_theta = 200;
phi = 0;
phi_dot = 0;
phi_ddot = 0;
%% Generalized coordinates 
q = [x;
    theta];
q_dot = [x_dot;
    theta_dot];
q_ddot = [x_ddot;
    theta_ddot];

%% Position vectors 
r2 = [x*cos(phi);
      L0 + x*sin(phi)];
r3 = [x*cos(phi) - cos(theta)*L2;
      L0 + x*sin(phi) - sin(theta)*L2] ;
r4 = [L1*cos(phi);
      L0 + L1*sin(phi)];

%% Velocity vectors
r_dot_2 = jacobian(r2,q)*q_dot; %+ jacobian(r2, phi)*phi_dot; %I am not sure if this is how to add the r,t part 
r_dot_3 = jacobian(r3,q)*q_dot; %+ jacobian(r3, phi)*phi_dot;
r_dot_4 = jacobian(r4,q)*q_dot; %+ jacobian(r4, phi)*phi_dot;

%% Kinetic energy 
T = (m2*x_dot^2)/2 + (m3*x_dot^2)/2 + (L1^2*m1*phi_dot^2)/6 + (L2^2*m3*theta_dot^2)/2 ... ;
    + (m2*phi_dot^2*x^2)/2 + (m3*phi_dot^2*x^2)/2 + L2*m3*theta_dot*x_dot*sin(theta - phi)... ;
    - L2*m3*phi_dot*theta_dot*x*cos(theta - phi);
%% Potential energy 
V = 1/2 * Lx^2 * kx + 1/2 * k_theta * theta^2 + 1/2* k_theta * theta_ref^2 + 1/2* k_theta * phi^2 ...;
    +1/2 * kx * x^2 + 1/2 * L0 * g * m0 + L0 * g *m1 + L0 * g *m2 + L0*g*m3 - Lx*kx *x - k_theta * theta * theta_ref ...;
    - k_theta * theta * phi + k_theta * theta_ref * phi - L2 * g * m3 * sin(theta) + 1/2 * L1 * g * m1 *sin(phi)...;
    + g * m2 * x * sin(phi) + g * m3 * x * sin(phi);
%% Generalized non-conservative force
Q_nc = [FA - dx * x_dot - dW* theta_dot * sin(theta - phi);
        d_theta * phi_dot - d_theta * theta_dot - L2 * dW * theta_dot]; 
%% Lagrange Equations of Motion
% d/dt(T,qdot) 
dTdq_dot = jacobian(T,q_dot);
First_Term1 = jacobian(transpose(dTdq_dot),q)*q_dot + jacobian(transpose(dTdq_dot),q_dot)*q_ddot;
              %+jacobian(transpose(dTdq_dot),phi)*phi_dot %+ jacobian(transpose(dTdq_dot),phi_dot)*phi_ddot;
First_Term = transpose(First_Term1); 

% T,q
Second_Term = jacobian(T,q);
Second_Term = simplify(Second_Term);

% V,q
Third_Term = jacobian(V,q);
Third_Term = simplify(Third_Term);

% Equation of Motion
EoM = First_Term - Second_Term + Third_Term - transpose(Q_nc);
EoM = simplify(EoM);
EoM = transpose(EoM);
%EoM == 0;


% Run Simulink model
simOut = sim('Simulink_deliverable_1');

% Time vector
t = simOut.tout;          % [0 110] s

% States q = [x; theta], stored as 1x2xN
q_raw = simOut.q;                 % 1x2x110001
q_mat = squeeze(q_raw).';         % -> 110001x2

x     = q_mat(:,1);               % x(t): hoist block position
theta = q_mat(:,2);               % θ(t): absolute chain angle


%% Equilibrium points 
dVdq = transpose(jacobian(V,q));
dVdq = simplify(dVdq);
q_01 = [12.5; -1.577];
q_02 = [12.5; 1.564];
q_01 = double(q_01);
q_02 = double(q_02);

%% Finding stability
% Find hessian and sub in equilibrium point
dVdq2 = hessian(V,q);
dVdq2 = simplify(dVdq2);
dVdq2_1 = subs(dVdq2,q,q_01); % sub in q0
dVdq2_2 = subs(dVdq2,q,q_02); % sub in q0

% Determine if positive definite
eig_val1 = eig(dVdq2_1);
eig_val2 = eig(dVdq2_2);
isposdef = all(eig_val1 > 0) % if =1 then stable.
isposdef = all(eig_val2 > 0) % if =1 then stable.

%% Linearization
M_0 = hessian(T,q_dot);
M_0 = subs(M_0,q,q_02); %sub in q = q_0
M_0 = subs(M_0,q_dot,[0;0]) % sub in q_dot = 0

D_0 = jacobian(-Q_nc,q_dot);
D_0 = subs(D_0,q,q_02);
D_0 = subs(D_0,q_dot,[0;0])

K_0 = hessian(V,q);
K_0 = subs(K_0,q,q_02);
K_0 = subs(K_0,q_dot,[0;0])

K_0_Q = jacobian(-Q_nc,q);
K_0_Q = subs(K_0_Q,q,q_02);
K_0_Q = subs(K_0_Q,q_dot,[0;0])

Q = subs(Q_nc,q,q_02);
Q = subs(Q,q_dot,[0;0])

EoM_linearized = M_0*q_ddot + D_0*q_dot + (K_0+K_0_Q)*(q-q_02) == Q

%% parameters
g_val        = 9.81;
L2_val       = 10;
m2_val       = 2000;
m3_val       = 500;
dx_val       = 200;
dW_val       = 2000;
d_theta_val  = 10;
kx_val       = 7000;
k_theta_val  = 200;
Lx_val       = 12.5;
phi_val      = 0;   % spd(t) = 0


M0 = double(subs(M_0, ...
    {g, L2, m2, m3, dx, dW, d_theta, kx, k_theta, Lx, phi}, ...
    {g_val, L2_val, m2_val, m3_val, dx_val, dW_val, d_theta_val, ...
     kx_val, k_theta_val, Lx_val, phi_val}));

D0 = double(subs(D_0, ...
    {g, L2, m2, m3, dx, dW, d_theta, kx, k_theta, Lx, phi}, ...
    {g_val, L2_val, m2_val, m3_val, dx_val, dW_val, d_theta_val, kx_val, k_theta_val, Lx_val, phi_val}));

K0 = double(subs(K_0, ...
    {g, L2, m2, m3, dx, dW, d_theta, kx, k_theta, Lx, phi}, ...
    {g_val, L2_val, m2_val, m3_val, dx_val, dW_val, d_theta_val, kx_val, k_theta_val, Lx_val, phi_val}));

KQ0 = double(subs(K_0_Q, ...
    {g, L2, m2, m3, dx, dW, d_theta, kx, k_theta, Lx, phi}, ...
    {g_val, L2_val, m2_val, m3_val, dx_val, dW_val, d_theta_val, kx_val, k_theta_val, Lx_val, phi_val}));

% Total stiffness matrix
Ktot = K0 + KQ0;

A = [zeros(2) eye(2);
     -M0\Ktot   -M0\D0];


% initial conditions
q0      = q_02;              
q_init  = [11.5; pi/3];      
q1_0    = q_init - q0;       
q1dot_0 = [0; 0];            

%state vector
z0 = [q1_0; q1dot_0];

%ode45
tspan = [0 110];
lin_ode = @(t,z) A*z;       

[t_lin, z_lin] = ode45(lin_ode, tspan, z0);

%something weird
q1_lin    = z_lin(:,1:2);
x_lin     = q1_lin(:,1) + q0(1);
theta_lin = q1_lin(:,2) + q0(2);

%% SIMULINK 
simOut = sim('Simulink_deliverable_1');

t_nl  = simOut.tout;
q_raw = simOut.q;         
q_mat = squeeze(q_raw).';  

x_nl     = q_mat(:,1);
theta_nl = q_mat(:,2);

%% Plot

%Plot x(t) nonlinear vs linearized
figure
plot(t_nl, x_nl, 'k', 'LineWidth', 1.2); hold on
plot(t_lin, x_lin, 'r', 'LineWidth', 1.2);
grid on
xlim([0 110]);
ylim([11.42 13.08]);       
xlabel('Time t [s]');
ylabel('Hoist position x(t) [m]');
title('Nonlinear vs linearized response of x(t)');
legend('Nonlinear system','Linearized system','Location','best');

%Plot theta(t) nonlinear vs linearized 
figure
plot(t_nl, theta_nl, 'k', 'LineWidth', 1.2); hold on
plot(t_lin, theta_lin, 'r', 'LineWidth', 1.2);
grid on
xlim([0 110]);
ylim([0.95 2.05]);        
xlabel('Time t [s]');
ylabel('\theta(t) [rad]');
title('Nonlinear vs linearized response of \theta(t)');
legend('Nonlinear system','Linearized system','Location','best');

