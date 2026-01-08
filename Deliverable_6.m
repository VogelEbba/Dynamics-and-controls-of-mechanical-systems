clc; delete all; clear all;

s = tf('s');
dw = 100; %For question a at least
G = (5.1e-3*s^2 + (0.2 + 0.049*dw)*s + 250)/ ... 
    (5e5*s^4 + (5e4 + 100*dw)*s^3 + 6.2e5*s^2 + 4.9e4*s);

%% Question a 
%C1 = K*((alfa/wc)*s +1)/((1/(alfa*wc))*s + 1);