clc; clear all; close all; 
syms theta theta_dot theta_ddot 
syms x x_dot x_ddot
syms phi phi_dot phi_ddot
syms L0 L1 L2 m0 m1 m2 m3 
syms g dx dW d_phi d_theta k_theta kx x0 theta0 q0 FA FW theta_ref Lx

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
r_dot_2 = jacobian(r2,q)*q_dot + jacobian(r2, phi)*phi_dot; %I am not sure if this is how to add the r,t part 
r_dot_3 = jacobian(r3,q)*q_dot + jacobian(r3, phi)*phi_dot;
r_dot_4 = jacobian(r4,q)*q_dot + jacobian(r4, phi)*phi_dot;

%% Kinetic energy 
T_crane_beam = 0;

%T_load = 1/2 * m3 * ( ... ;
%    (x_dot * cos(phi) + theta_dot * sin(theta) * L2 - x * sin(phi) * phi_dot)^2 + ... ;
%    (x_dot * sin(phi) - theta_dot * cos(theta) * L2 + x * cos(phi) * phi_dot)^2);

%T_jib = 1/2* m1 *(1/4* L1^2 * phi_dot^2 * cos(phi)^2 + 1/4 * L1^2 * phi_dot^2 * sin(phi)^2)...
%        + 1/24 * m1 * L1^2 * phi_dot^2; 

%T_hoist_block  = 1/2 *m2 *((x_dot*cos(phi) - x*sin(phi)*phi_dot)^2 ... ;
%                 + (x_dot*sin(phi) + x*cos(phi)*phi_dot)^2 ); 

%T = T_crane_beam + T_load + T_jib + T_hoist_block;
%T = simplify(T);
T = (m2*x_dot^2)/2 + (m3*x_dot^2)/2 + (L1^2*m1*phi_dot^2)/6 + (L2^2*m3*theta_dot^2)/2 ... ;
    + (m2*phi_dot^2*x^2)/2 + (m3*phi_dot^2*x^2)/2 + L2*m3*theta_dot*x_dot*sin(theta - phi)... ;
    - L2*m3*phi_dot*theta_dot*x*cos(theta - phi);
%% Potential energy 
%V_load = m3 * g * (L0 + x * sin(phi) - L2 * sin(theta)) + 1/2 * k_theta * (theta_ref - theta + phi)^2;
%V_jib = m1 * g * 0.5 * L1 * sin(phi + L0);
%V_hoist_block = m2 * g * (x * sin(phi) + L0) + (1/2) * kx * (x - x0)^2;
%V_crane_beam = m0 * g * (1/2) * L0;
%V = V_load + V_jib + V_hoist_block + V_crane_beam;
V = 1/2 * Lx^2 * kx + 1/2 * k_theta * theta^2 + 1/2* k_theta * theta_ref^2 + 1/2* k_theta * phi^2 ...;
    +1/2 * kx * x^2 + 1/2 * L0 * g * m0 + L0 * g *m1 + L0 * g *m2 + L0*g*m3 - Lx*kx *x - k_theta * theta * theta_ref ...;
    - k_theta * theta * phi + k_theta * theta_ref * phi - L2 * g * m3 * sin(theta) + 1/2 * L1 * g * m1 *sin(phi)...;
    + g * m2 * x * sin(phi) + g * m3 * x * sin(phi);
%% Generalized non-conservative force
%Force part
%Force FA
%r_FA= [x * cos(phi);
%      L0 + x * sin(phi)];
%FA_vec = [cos(phi)*FA;
%          sin(phi)*FA];
%Force dx 
%r_dx = r2;
%F_dx = [-cos(phi)*x_dot*dx;
%        -sin(phi)*x_dot*dx];

%Moment part 
%Moment d theta
%angle_d_theta = [0;  %Used angle instead of theta because theta is already used as a generalized coordinate. 
%         0;
%         theta]; 
%M_d_theta= [0
%        0
%        -theta_dot* d_theta];
%moment voor d phi 
%angle_d_spd = [0;   
%               0;
%M_d_spd = [0 
%           0
%           -phi_dot * d_phi]; %Deze d_phi can ook naar d_spd worden verandert. 
%Moment for dw
%angle_dw = [0                 %Deze is nog onduidelijk of hij klopt
%            0
%            theta];
%M_dw = [0
%        0
%        -dW*theta_dot];


%Q_nc = transpose(jacobian(r_FA,q))*FA_vec + transpose(jacobian(r_dx,q))*F_dx ...;
%        + transpose(jacobian(angle_d_theta,q))*M_d_theta...;
%        + transpose(jacobian(angle_d_spd,q))*M_d_spd...;
%        + transpose(jacobian(angle_dw,q))*M_dw;

Q_nc = [FA - dx * x_dot - dW* theta_dot * sin(theta - phi);
        d_theta * phi_dot - d_theta * theta_dot - L2 * dW * theta_dot]; 
%% Lagrange Equations of Motion
% d/dt(T,qdot) 
dTdq_dot = jacobian(T,q_dot);
First_Term1 = jacobian(transpose(dTdq_dot),q)*q_dot + jacobian(transpose(dTdq_dot),q_dot)*q_ddot...;
              +jacobian(transpose(dTdq_dot),phi)*phi_dot + jacobian(transpose(dTdq_dot),phi_dot)*phi_ddot;
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

%% Plot 1: x(t)
figure
plot(t, x)
grid on
xlim([0 110])
ylim([11.15 12.26])
xlabel('Time t [s]')
ylabel('Hoist position x(t) [m]')
title('Time trajectory of the hoist block position x(t)')
legend('x(t)')

%% Plot 2: θ(t)
figure
plot(t, theta)
grid on
xlim([0 110])
ylim([0.95 2.05])
xlabel('Time t [s]')
ylabel('\theta(t) [rad]')
title('Time trajectory of the absolute chain angle \theta(t)')
legend('\theta(t)')
