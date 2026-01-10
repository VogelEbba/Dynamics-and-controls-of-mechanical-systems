clc; delete all; clear all;

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