clear all
clc
close all
s = tf('s');
%% Computing the poles of the transfer function
Gl = tf([5.1e-3 100 250 ],[5e5 2.5e5 6.4e5 4.9e4 0]); 
H = tf([5.1e-3 100 250 0 ],[5e5 2.5e5 6.4e5 4.9e4 0]);


Hmin = minreal(H);      % this cancels common factors
poles_Hmin = pole(Hmin);

%% Generating impuls response function
t = 0:0.01:91; 
[y, tOut] = impulse(Hmin,t);
figure(1);
plot(tOut, y);
grid on
xlim([0 91]);
xlabel('time [s]');
ylabel('Amplitude');

%% Stability of the closed loop transfer function 
P1 = 80;
P2 = 400;
P3 = 700;

L1 = tf([P1*5.1e-3 P1*100 P1*250], [5e5, 2.5e5, 6.4e5, 4.9e4, 0]);
L2 = tf([P2*5.1e-3 P2*100 P2*250], [5e5, 2.5e5, 6.4e5, 4.9e4, 0]);
L3 = tf([P3*5.1e-3 P3*100 P3*250], [5e5, 2.5e5, 6.4e5, 4.9e4, 0]);

T1_feedback = feedback(L1, 1);
T2_feedback = feedback(L2, 1);
T3_feedback = feedback(L3, 1);

% T1 = tf([P1*5.1e-3 P1*100 P1*250],[5e5 2.5e5 (6.4e5+P1*5.1e-3) (4.9e4 + P1* 100) P1*250]);
% T2 = tf([P2*5.1e-3 P2*100 P2*250],[5e5 2.5e5 (6.4e5+P2*5.1e-3) (4.9e4 + P2* 100) P2*250]);
% T3 = tf([P3*5.1e-3 P3*100 P3*250],[5e5 2.5e5 (6.4e5+P3*5.1e-3) (4.9e4 + P3* 100) P3*250]);

poles_T1 = pole(T1_feedback);
poles_T2 = pole(T2_feedback);
poles_T3 = pole(T3_feedback);




%% Step response function 
t2 = 0 : 0.01 : 500;
[Amplitude_step1, tOut_step1] = step(T1_feedback,t2);
[Amplitude_step2, tOut_step2] = step(T2_feedback,t2);
[Amplitude_step3, tOut_step3] = step(T3_feedback,t2);

figure(2);
plot(tOut_step1, Amplitude_step1);
hold on 
plot(tOut_step2, Amplitude_step2);
hold on 
plot(tOut_step3, Amplitude_step3);
legend("P1", "P2","P3")
title("Step response of T(s)")
grid on
xlim([0 240]);
xlabel('time [s]');
ylabel('Amplitude');

%% poles and zeros of the closed-loop transfer functions
figure(3);
pzmap(T1_feedback); 
hold on;

pzmap(T2_feedback);
hold on;

pzmap(T3_feedback);
grid on;

xlim([-0.24 -5.5e-3]);
ylim([-1.2 1.2]);
title('Pole and zeros of closed-loop transfer functions');
legend('P1','P2','P3','Location','best');

% Zoomed in verion of ths plot
figure(4);
pzmap(T1_feedback); 
hold on;

pzmap(T2_feedback);
hold on;

pzmap(T3_feedback);
grid on;

xlim([-0.05 -0.031]);
ylim([-0.74 0.74]);
title('Zoomed in version');
legend('P1','P2','P3','Location','best');

%% Step response of the closed-loop system from r to y
C_PD = tf([390 360], 1);
T_PD = feedback(C_PD*Gl,1);

t_pd = 0:0.01:79;

[y_pd, t_pd_out] =step(T_PD,t_pd);     
[y_P2, t_P2] =step(T2_feedback,t_pd);    

figure(5);
plot(t_P2, y_P2); 
hold on;

plot(t_pd_out, y_pd);
grid on;

xlim([0 79]);
xlabel('time (s)');
ylabel('amplitude');
title('Step response comparison');
legend('P2 = 400','C(s) = 390s + 360','Location','best');

%% location of poles and zeros for both cases
figure(6);
pzmap(T2_feedback); 
hold on;

pzmap(T_PD);
grid on;

xlim([-0.22 -0.025]);
ylim([-1.2 1.2]);
title('location poles and zeros');
legend('P2 = 400','PD controller','Location','best');

poles_P2  = pole(T2_feedback);
poles_PD  = pole(T_PD);

%% Final value
y1_step_FV = evalfr(T1_feedback,0);
finalvalueTPD = evalfr(T_PD,0);
S_T1 = stepinfo(T1_feedback, 'SettlingTimeThreshold', 0.01)
S_T2 = stepinfo(T2_feedback, 'SettlingTimeThreshold', 0.01)
S_T3 = stepinfo(T3_feedback, 'SettlingTimeThreshold', 0.01)
S_PD = stepinfo(T_PD, 'SettlingTimeThreshold', 0.01);
