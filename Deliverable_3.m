clc; clear all; close all;

syms theta theta_dot theta_ddot 
syms x x_dot x_ddot t s

%% Parameters 
m1 = 10000;
m2 = 2000;
m3 = 500;
L2 = 10;
kx = 0;
k_theta = 200;
g = 9.81;
dW = 2000;
d_theta = 10;
dx = 200;

%% Generalized coordinates 
q = [x;
    theta];
q_dot = [x_dot;
    theta_dot];
q_ddot = [x_ddot;
    theta_ddot];

%% Stable equilibrium 
q_0 = [12.5;
       1.564];

%% M_0 and K_0 (From the solutions of deliverable 2)
M_0 = [m2+m3, -L2*m3*sin(-1.564);
       -L2*m3*sin(-1.564), L2^2*m3];
K_0 = [kx, 0;
       0, k_theta + L2*g*m3];
D_0 = [dx, -dW*sin(-1.564);
       0, d_theta + L2*dW];

%% Eigenmodes and eigenfrequencies

[U,lambda] = eig(M_0\K_0);     % Solve for eigenmatrix (U) and eigenvalues matrix (lambda)
u1 = U(:,1);                    % 1st eigenmode
u2 = U(:,2);                    % 2nd eigenmode


u1_norm = u1/u1(1);
u2_norm = u2/u2(1);

w = sqrt(diag(lambda));       % Solve for eigenfrequencies (rad/s)
w1 = w(1);                    % 1st eigenfrequency
w2 = w(2);                   % 2nd eigenfrequency

%% SIMULINK (nonlinear system)
% simOut = sim('Simulink_deliverable_1');
% 
% t_nl  = simOut.tout;
% q_raw = simOut.q;         
% q_mat = squeeze(q_raw).';  
% 
% x_nl     = q_mat(:,1);
% theta_nl = q_mat(:,2);
% 
% %% Plot
% 
% %xt nonlinear vs linearized
% figure
% plot(t_nl, x_nl, 'k', 'LineWidth', 1); hold on      % nonlinear (Simulink)
% grid on
% xlim([0 130]);
% ylim([11.5 22.5]);       
% xlabel('Time t [s]');
% ylabel('Hoist position x(t) [m]');
% title('Nonlinear response of x(t)');
% 
% %theta t nonlinear vs linearized 
% figure
% plot(t_nl, theta_nl, 'k', 'LineWidth', 1.2); hold on  % nonlinear (Simulink)
% grid on
% xlim([0 130]);
% ylim([1.38 1.821]);        
% xlabel('Time t [s]');
% ylabel('\theta(t) [rad]');
% title('Nonlinear response of \theta(t)');
 
%% Computing damped transfer functions
B = [1;
     0];

A = [zeros(2) eye(2);
     -M_0\K_0   -M_0\D_0];

Bss = [zeros(2,1);
       M_0\B];

C = [1 0 0 0;
     0 1 0 0];

D = [0;0];
I = eye(4)

sys = ss(A,Bss,C,D);
G = tf(sys) % gives G1(s) and G2(s)
G1 = G(1);
G2 = G(2);
figure;
bode(G)  

%% Gain of the transferfunctions
[mag, phase] = bode(G1, w2);
[mag,~] = bode(G1, w2);
gain1_option3 = mag(:)

[mag, phase] = bode(G2, w2);
[mag,~] = bode(G2, w2);
gain2_option3 = mag(:)

%% Question f

Gl = G1 + (L2 * G2);

figure(1);
bode (Gl)