%% Parameters 
m1 = 10000;
m3 = 500;
L2 = 10;
kx = 0;
k_theta = 200;
g = 9.81;

%% M_0 and K_0 (From the solutions of deliverable 2)
M_0 = [m1+m3, -L2*m3*sin(-1.564);
       -L2*m3*sin(-1.564), L2^2*m3];
K_0 = [kx, 0;
       0, k_theta + L2*g*m3];

%% Eigenmodes and eigenfrequencies

[U,lambda] = eig(M_0\K_0);     % Solve for eigenmatrix (U) and eigenvalues matrix (lambda)
u1 = U(:,1)                    % 1st eigenmode
u2 = U(:,2)                    % 2nd eigenmode

w = sqrt(diag(lambda));       % Solve for eigenfrequencies (rad/s)
w1 = w(1)                     % 1st eigenfrequency
w2 = w(2)                     % 2nd eigenfrequency
 