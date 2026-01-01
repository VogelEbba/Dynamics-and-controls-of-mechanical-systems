clc;
clear all;
close all;

s = tf('s');
G = tf([5.1e-3 100 250],[5e5 2.5e5 6.4e5 4.9e4 0]);

%% a
C1 = (770*s + 350)/(4.3*s + 1);
C2 = (770*s + 210)/(2.5*s + 1);
C3 = (770*s + 140)/(1.7*s + 1);

L1 = C1 * G;
L2 = C2 * G;
L3 = C3 * G;

figure(1);
nyquist(L1); hold on
nyquist(L2); hold on
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

MM1 = 1 / norm(S1, inf); %Twijfels of dit goed is
MM2 = 1 / norm(S2, inf);
MM3 = 1 / norm(S3, inf);

% Andere manier: Geeft hetzelfde antwoord
w = linspace(0,1000,1e6);
[modulusmargin1,index1] = min(abs(freqresp(1 + L1,w)));
[modulusmargin2,index2] = min(abs(freqresp(1 + L2,w)));
[modulusmargin3,index3] = min(abs(freqresp(1 + L3,w)));
wc1 = w(index1);
wc2 = w(index2);
wc3 = w(index3);

%% d step response
t = 0 : 0.01 : 150;
[Amplitude_step1, tOut_step1] = step(T1 ,t);
[Amplitude_step2, tOut_step2] = step(T2 ,t);
[Amplitude_step3, tOut_step3] = step(T3,t);

figure(2);
plot(tOut_step1, Amplitude_step1); hold on
plot(tOut_step2, Amplitude_step2); hold on
plot(tOut_step3, Amplitude_step3);
legend("C1", "C2","C3")
title("Step response of T(s)")
grid on
ylim([-4 4])
xlim([0 100]);
xlabel('time [s]');
ylabel('Amplitude');

%info1 = stepinfo(T1, 'SettlingTimeThreshold', 0.01);
info2 = stepinfo(T2, 'SettlingTimeThreshold', 0.01);
info3 = stepinfo(T3, 'SettlingTimeThreshold', 0.01);

%% f
G_damaged = tf([5.1e-3 10 250],[5e5 7e4 6.2e5 4.9e4 0]);

L1_d = C1 * G_damaged;
L2_d = C2 * G_damaged;
L3_d = C3 * G_damaged;

T1_d = feedback(L1_d, 1);
T2_d = feedback(L2_d, 1);
T3_d = feedback(L3_d, 1);

t = 0 : 0.01 : 150;

figure;
step(T1_d, T2_d, T3_d, t)
grid on
xlim([0 100])
ylim([-4 4])
xlabel('time [s]');
ylabel('Amplitude');
title('Step response with damaged damping (90% decrease)');
legend("C1 damaged", "C2 damaged", "C3 damaged", 'location', 'best')

info1_d = stepinfo(T1_d,'SettlingTimeThreshold',0.01);
info2_d = stepinfo(T2_d,'SettlingTimeThreshold',0.01);
info3_d = stepinfo(T3_d,'SettlingTimeThreshold',0.01);


%% g
info1 = stepinfo(T1,'SettlingTimeThreshold',0.01);
info2 = stepinfo(T2,'SettlingTimeThreshold',0.01);
info3 = stepinfo(T3,'SettlingTimeThreshold',0.01);

fprintf('--- Settling time original ---\n');
fprintf('C1: %.2f s, C2: %.2f s, C3: %.2f s\n', info1.SettlingTime, info2.SettlingTime, info3.SettlingTime);

fprintf('--- Settling time damaged ---\n');
fprintf('C1: %.2f s, C2: %.2f s, C3: %.2f s\n', info1_d.SettlingTime, info2_d.SettlingTime, info3_d.SettlingTime);

%% h
S1 = 1/(1 + L1);
S2 = 1/(1 + L2);
S3 = 1/(1 + L3);

%modulus margin
MM1 = 1 / norm(S1, inf);
MM2 = 1 / norm(S2, inf);
MM3 = 1 / norm(S3, inf);

fprintf('Modulus margins: MM1=%.4f, MM2=%.4f, MM3=%.4f\n', MM1, MM2, MM3);

figure;
bodemag(S1, S2, S3);
grid on
title('Sensitivity magnitude |S(j\omega)| (original d_W)');
legend('S1','S2','S3','Location','best');

%% i
figure;
bodemag(T1, T2, T3);
grid on
title('Complementary sensitivity magnitude');
legend('T1','T2','T3','Location','best');

%% j
figure;
bodemag(S2, T2);
grid on
title('Bode of S2(s) and T2(s)');
legend('S2','T2','Location','best');



