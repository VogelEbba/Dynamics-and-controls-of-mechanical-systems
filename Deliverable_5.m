clc; clear all; delete all 
s = tf('s');
G = tf([5.1e-3 100 250 ],[5e5 2.5e5 6.4e5 4.9e4 0]); 

%% a 
C1 = (770*s + 350)/(4.3*s + 1);
C2 = (770*s + 210)/(2.5*s + 1);
C3 = (770*s + 140)/(1.7*s + 1);

L1 = C1 * G;
L2 = C2 * G;
L3 = C3 * G;

figure;
nyquist(L1);
hold on 
nyquist(L2);
hold on 
nyquist(L3);
xlim([-2 0.2])
ylim([-1 1])
legend('L1','L2','L3')
title('Nyquist plot for L1, L2 and L3')

%% b
T1 = feedback(L1, 1);
T2 = feedback(L2, 1);
T3 = feedback(L3, 1);

L1_poles = pole(L1);
L2_poles = pole(L2);
L3_poles = pole(L3);

%% c 
S1 = 1/(1 + L1);
S2 = 1/(1 + L2);
S3 = 1/(1 + L3); 

[Gm_L1, Pm_L1, Wcg_L1, Wcp_L1] = margin(L1);
[Gm_L2, Pm_L2, Wcg_L2, Wcp_L2] = margin(L2);
[Gm_L3, Pm_L3, Wcg_L3, Wcp_L3] = margin(L3);

MM1 = 1 / norm(S1, inf);
MM2 = 1 / norm(S2, inf);
MM3 = 1 / norm(S3, inf);