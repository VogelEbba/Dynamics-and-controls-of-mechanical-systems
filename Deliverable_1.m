clc; clear all; close all; 
syms theta theta_dot theta_ddot 
syms x x_dot x_ddot
syms s_pd s_pd_dot s_pd_ddot
syms L0 L1 L2 m0 m1 m2 m3 
syms g dx dW d_phi d_theta k_theta kx x0 theta0 q0 FA FW

%% Generalized coordinates 
q = [x;
    theta];
q_dot = [x_dot;
    theta_dot];
q_ddot = [x_ddot;
    theta_ddot];

%% Position vectors 
r2 = [x*cos(s_pd);
      L0 + x*sin(s_pd)];
r3 = [x*cos(s_pd) - cos(theta)*L2;
      L0 + x*sin(s_pd) - sin(theta)*L2] ;
r4 = [L1*cos(s_pd);
      L0 + L1*sin(s_pd)];

%% Velocity vectors
r_dot_2 = jacobian(r2,q)*q_dot + jacobian(r2, s_pd)*s_pd_dot; %I am not sure if this is how to add the r,t part 
r_dot_3 = jacobian(r3,q)*q_dot + jacobian(r3, s_pd)*s_pd_dot;
r_dot_4 = jacobian(r4,q)*q_dot + jacobian(r4, s_pd)*s_pd_dot;

%% Kinetic energy 
T_crane_beam = 0;
T_load = 1/2 * m3 * ( ... ;
    (x_dot * cos(s_pd) - theta_dot * sin(theta) * L2 - x * sin(s_pd) * s_pd_dot)^2 + ... ;
    (x_dot * sin(s_pd) + theta_dot * cos(theta) * L2 + x * cos(s_pd) * s_pd_dot)^2 ... ;
);
T_jib = 1/24 * m1 * L1 * s_pd_dot^2; 
T_hoist_block  = 1/2 *m2 *( ... ;
    (x_dot*cos(s_pd) - x*sin(s_pd)*s_pd_dot)^2 + ... ;
    (x_dot*sin(s_pd) + x*cos(s_pd)*s_pd_dot)^2 ); 

T = T_crane_beam + T_load + T_jib + T_hoist_block;
T = simplify(T);

%% Potential energy 
V_load = m3 * g * (L0 + x * sin(s_pd) - L2 * sin(theta)) + (1/2) * k_theta * (theta - theta0)^2;
V_jib = m1 * g * (1/2) * L1 * sin(s_pd + L0);
V_hoist_block = m2 * g * (x * sin(s_pd) + L0) + (1/2) * kx * (x - x0)^2;
V_crane_beam = m0 * g * (1/2) * L0;
V = V_load + V_jib + V_hoist_block + V_crane_beam;

%% Generalized non-conservative force
%Force part
%Force FA
r_FA= [x * cos(s_pd);
      L0 + x * sin(s_pd)];
FA_vec = [cos(s_pd)*FA;
          sin(s_pd)*FA];
%Force dx 
r_dx = r2;
F_dx = [-cos(s_pd)*x_dot*dx;
        -sin(s_pd)*x_dot*dx];

%Moment part 
%Moment d theta
angle_d_theta = [0;  %Used angle instead of theta because theta is already used as a generalized coordinate. 
         0;
         theta]; 
M_d_theta= [0
        0
        -theta_dot* d_theta];
%moment voor d phi 
angle_d_spd = [0;   
               0;
               s_pd]; 
M_d_spd = [0 
           0
           -s_pd_dot * d_phi]; %Deze d_phi can ook naar d_spd worden verandert. 
%Moment for dw
angle_dw = [0                 %Deze is nog onduidelijk of hij klopt
            0
            theta];
M_dw = [0
        0
        -dW*theta_dot];


Q_nc = transpose(jacobian(r_FA,q))*FA_vec + transpose(jacobian(r_dx,q))*F_dx ...;
        + transpose(jacobian(angle_d_theta,q))*M_d_theta...;
        + transpose(jacobian(angle_d_spd,q))*M_d_spd...;
        + transpose(jacobian(angle_dw,q))*M_dw;

%% Lagrange Equations of Motion
% d/dt(T,qdot) 
dTdq_dot = jacobian(T,q_dot);
First_Term = jacobian(transpose(dTdq_dot),q)*q_dot + jacobian(transpose(dTdq_dot),q_dot)*q_ddot;
First_Term = transpose(First_Term); 

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
EoM == 0;


%out = sim("Simulink_deliverable_1.slx");

figure(1)
plot(out.tout, out.q(:,1),'k');
grid on
xlabel('Time [s]');
ylabel('X [m]');


figure(2)
plot(out.tout, out.q(:,2),'k');
grid on
xlabel('Time [s]');
ylabel('Theta [rad]');
