clc; delete all; 
close all;

s = tf('s');
dwa = 100; 
dwb = 10;
Ga = (5.1e-3*s^2 + (0.2 + 0.049*dwa)*s + 250)/ ... 
    (5e5*s^4 + (5e4 + 100*dwa)*s^3 + 6.2e5*s^2 + 4.9e4*s);
Gb = (5.1e-3*s^2 + (0.2 + 0.049*dwb)*s + 250)/ ... 
    (5e5*s^4 + (5e4 + 100*dwb)*s^3 + 6.2e5*s^2 + 4.9e4*s);
figure(1);
bode(Ga)
figure(2);
margin(Ga)
title('Margin G');
%% Question a 
fc = 0.025;    %Hz 
wc = fc*2*pi; %rad/s
alpha = 1.2; % This is a guess

C0 = ((alpha/wc)*s +1)/((1/(alpha*wc))*s + 1); %C1 without K
[mag,phase,wout] = bode(C0*Ga, wc);
K = 1/mag;

C1 = K * C0;
L1 = C1 * Ga;
L1 = minreal(L1);
S1 = 1/(1 + L1);
T1 = feedback(L1, 1);

MM1 = norm(S1, inf);
MM1_dB = 20*log10(MM1);

figure(3);
nyquist(L1)
ylim([-5 5])
xlim([-2 1])
L1_poles = pole(L1); %L1 is stabiel als je er vanuit gaat dat een pole die 0 is kan negeren, ik heb dit op de discussie pagina gevraagd

%% Question b 
[Gm_L1, Pm_L1, Wcg_L1, Wcp_L1] = margin(L1);

figure(4);
margin(L1)
title('Margin L1'); 

%% Question c 
t = 0 : 0.1 : 100;
[y, t_out] = step(T1, t);

r = ones(size(t));
e = r - y;

figure(5);
plot(t_out, y);
title('Step response of the closed loop system');
xlim([0 100]);
xlabel('time (s)');
ylabel('xl (m)');

figure(6);
plot(t_out, e);
title('Step response of the error');
xlim([0 100]);
xlabel('time (s)');
ylabel('e (m)');

figure(7);
plot(t_out, y, 'Color', 'b');
hold on
plot(t_out, e, 'Color','r');
title('The output of x_l and e combined');
legend('xl','e');
xlabel('time (s)');
ylabel('e (m) and xl(m)');


%De code voor dit stuk is grotendeels met chat geschreven. 
tol = 0.025;
dt = t_out(2) - t_out(1);
window = round(5/dt);

insideBand = abs(e) < tol;
settlingTime = NaN;

for k = 1:length(insideBand) - window
    if all(insideBand(k:k+window))
        settlingTime = t_out(k);
        break
    end
end

settlingTime;

%% Question d

L2 = C1 * Gb;
L2 = minreal(L2);
S2 = 1/(1 + L2);
T2 = feedback(L2, 1);

MM2 = norm(S2, inf);
MM2_dB = 20*log10(MM2);

n = 1000; 
theta = linspace(0, 2*pi, n);
x = cos(theta);
y = sin(theta);

figure(8);
nyquist(L1);
hold on
nyquist(L2);
hold on
plot(x, y, 'Color', 'black');
legend('L1','L2','Unit circle');
ylim([-2, 2]);
xlim([-2.5, 0.5]);

L2_poles = pole(L2); 

[Gm_L2, Pm_L2, Wcg_L2, Wcp_L2] = margin(L2);

figure(9);
margin(L2)
title('Margin L2'); 

%% Question e
s = tf('s');

%given
fc2 = 0.06;            
wc2 = 2*pi*fc2;        

dW = 100;
G2 = (5.1e-3*s^2 + (0.2 + 0.049*dW)*s + 250) / ...
     (5e5*s^4 + (5e4 + 100*dW)*s^3 + 6.2e5*s^2 + 4.9e4*s);

alpha = 3;

wz = 1.11;    
wp = 1.30;     

bz = 0.02;      
bp = 0.3;     


Clead = ((alpha/wc2)*s + 1) / (s/(alpha*wc2) + 1);
Cnotch = (wp/wz)^2 * (s^2 + 2*bz*wz*s + wz^2) / (s^2 + 2*bp*wp*s + wp^2);

C2_woK = minreal(Clead * Cnotch);

%|L(j*wc2)| = 1 
magL = squeeze(bode(C2_woK * G2, wc2));
K2   = 1/magL;

C2 = minreal(K2 * C2_woK);
L2e = minreal(C2 * G2);       

% Closed-loop and sensitivity
T2e = feedback(L2e, 1);
S2e = 1/(1 + L2e);

MM2    = norm(S2e, inf);
MM2_dB = 20*log10(MM2);


fprintf('K2 = %.4g\n', K2);
fprintf('MM = %.3f (%.2f dB)\n', MM2, MM2_dB);
fprintf('Closed-loop stable? %d\n', isstable(T2e));

figure(20); clf;
nyquist(L2e); 
grid on;
xlim([-2 1]); 
ylim([-5 5]);
title('Nyquist contour of L2(s) (Question e)');

figure(21); clf;
bodemag(L2e); 
grid on;
title('Magnitude plot of L2(s) (Question e)');
xline(wc2,'r--','\omega_c','LabelVerticalAlignment','bottom');
yline(1,'k--','|L| = 1');

figure(22); clf;
bodemag(S2e); 
grid on;
title('Magnitude plot of S2(s) (Question e)');
yline(6,'k--','6 dB');

%% Question f 
t = 0:0.1:100;

y1 = step(T1, t);   
y2 = step(T2e, t);   

% errors
e1 = 1 - y1;
e2 = 1 - y2;

figure; 
clf;
plot(t, y1); 
hold on;
plot(t, y2);
grid on; 
xlim([0 100]);
xlabel('time [s]'); 
ylabel(['x_l [m]']);
legend('C1','C2','Location','best');
title('x_l(t): C1 vs C2');

figure; 
clf;
plot(t, e1); 
hold on;
plot(t, e2);
yline(0.025); 
yline(-0.025);
grid on; 
xlim([0 100]);
xlabel('time [s]'); 
ylabel('e(t) [m]');
legend('e (C1)','e (C2)','\pm 0.025','Location','best');
title('e(t): C1 vs C2');

tol = 0.025;
dt = t(2) - t(1);
window = round(5/dt);

inside = abs(e2) < tol;
t_set = NaN;
for k = 1:length(inside)-window
    if all(inside(k:k+window))
        t_set = t(k);
        break
    end
end

fprintf('Settling time (C2), |e(t)|<0.025 for 5s: %.2f s\n', t_set);


%% question g

figure; 
clf;
bodemag(S1); 
hold on;
bodemag(S2e);  

grid on;
legend('S1 (C1)','S2 (C2)','Location','best');
title('Sensitivity magnitude: S1(s) vs S2(s)');

%% Question h
Ar    = 0.1;        
wr    = 0.01;       
tEnd  = 1500;      
dt    = 0.5;      
t     = 0:dt:tEnd;

r = Ar*sin(wr*t);  

y1 = lsim(T1,  r, t);
y2 = lsim(T2e, r, t);

e1 = r(:) - y1;
e2 = r(:) - y2;

figure; 
clf;
plot(t, e1); 
hold on;
plot(t, e2);
grid on; 
xlim([0 tEnd]);
xlabel('time [s]');
ylabel('error e(t) [m]');
title('Tracking error for sin reference');
legend('e1 (C1)','e2 (C2)','Location','best');

%% Question i


Ad = 100;
wd = 0.8;          
t  = 0:0.01:50;    
d  = Ad*sin(wd*t);

Gd1 = minreal(G2 * S1);     
Gd2 = minreal(G2 * S2e);    

y1d = lsim(Gd1, d, t);
y2d = lsim(Gd2, d, t);

e1 = -y1d;
e2 = -y2d;

figure; 
clf;
plot(t, e1); 
hold on;
plot(t, e2);
grid on; 
xlim([0 50]);
xlabel('time [s]');
ylabel('error e(t) [m]');
title('Error due to input disturbance d(t): e1 vs e2');
legend('e1 (C1)', 'e2 (C2)', 'Location', 'best');