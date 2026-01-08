clc; delete all; clear all;

s = tf('s');
dw = 100; %For question a at least
G = (5.1e-3*s^2 + (0.2 + 0.049*dw)*s + 250)/ ... 
    (5e5*s^4 + (5e4 + 100*dw)*s^3 + 6.2e5*s^2 + 4.9e4*s);

%% Question a 
fc = 0.025;    %Hz 
wc = fc*2*pi; %rad/s
alpha = 0.3; % This is a guess

C0 = ((alpha/wc)*s +1)/((1/(alpha*wc))*s + 1); %C1 without K
[mag,phase,wout] = bode(C0*G, wc);
K = 1/mag;

C1 = K * C0;
L1 = C1 * G;
S1 = 1/(1 + L1);
MM1 = norm(S1, inf)

%Stability needs to be checked