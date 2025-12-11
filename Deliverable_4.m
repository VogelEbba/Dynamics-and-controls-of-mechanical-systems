clear all
clc
close all

%% Computing the poles of the transfer function
Gl = tf([5.1e-3 100 250 ],[5e5 2.5e5 6.4e5 4.9e4 0]); 
H = tf([5.1e-3 100 250 0 ],[5e5 2.5e5 6.4e5 4.9e4 0]);


Hmin = minreal(H);      % this cancels common factors
poles_Hmin = pole(Hmin);

%% Generating impuls response function
t = 0:0.01:91; 
[y, tOut] = impulse(H,t);
figure;
plot(tOut, y);
grid on
xlim([0 91]);
xlabel('time [s]');
ylabel('Amplitude');

%% Stability of the closed loop transfer function 
P1 = 80;
P2 = 400;
P3 = 700;


T1 = tf([P1*5.1e-3 P1*100 P1*250],[5e5 2.5e5 (6.4e5+P1*5.1e-3) (4.9e4 + P1* 100) P1*250]);
T2 = tf([P2*5.1e-3 P2*100 P2*250],[5e5 2.5e5 (6.4e5+P2*5.1e-3) (4.9e4 + P2* 100) P2*250]);
T3 = tf([P3*5.1e-3 P3*100 P3*250],[5e5 2.5e5 (6.4e5+P3*5.1e-3) (4.9e4 + P3* 100) P3*250]);
poles_T = pole(T1);
Tmin = minreal(T1);
poles_Tmin = pole(Tmin);

%% Step response function 
t2 = 0 : 0.01 : 500;
[Amplitude_step1, tOut_step1] = step(T1,t2);
[Amplitude_step2, tOut_step2] = step(T2,t2);
[Amplitude_step3, tOut_step3] = step(T3,t2);

figure;
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
figure;
pzmap(T1); 
hold on;

pzmap(T2);
hold on;

pzmap(T3);
grid on;

xlim([-0.24 -5.5e-3]);
ylim([-1.2 1.2]);
title('Pole and zeros of closed-loop transfer functions');
legend('P1','P2','P3','Location','best');

% Zoomed in verion of ths plot
figure;
pzmap(T1); 
hold on;

pzmap(T2);
hold on;

pzmap(T3);
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
[y_P2, t_P2] =step(T2,t_pd);    

figure;
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
figure;
pzmap(T2); 
hold on;

pzmap(T_PD);
grid on;

xlim([-0.22 -0.025]);
ylim([-1.2 1.2]);
title('location poles and zeros');
legend('P2 = 400','PD controller','Location','best');

poles_P2  = pole(T2);
poles_PD  = pole(T_PD);
