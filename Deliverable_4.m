clear all
clc
delete all
syms s
%% Computing the poles of the transfer function
Gl = tf([5.1e-3 1000 250 ],[5e5 2.5e5 6.4e5 4.9e4 0]);
H = tf([5.1e-3 1000 250 0 ],[5e5 2.5e5 6.4e5 4.9e4 0]);


Hmin = minreal(H);      % this cancels common factors
poles_Hmin = pole(Hmin)

%% Generating impuls response function
t = 0:0.01:91; 
[y, tOut] = impulse(H,t);
figure;
plot(tOut, y);
grid on
xlim([0 91]);
xlabel('time [s]');
ylabel('Amplitude');

