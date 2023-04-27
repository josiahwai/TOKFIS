% define kp, ki, kd coefficients for the PID controller
% the coefficients correspond with tok.ccnames

function [K, Kp, Ki, Kd] = nstxu_coilcontrol_pid(tok)

% a simple method is for the controller weights to be proportional to the
% circuit inductance:
L = diag(tok.mcc);

kp = L * 120;    
ki = kp * 20;
kd = kp * 0;     % set the derivative term to zero
tauroll = 1e-2;

s = tf('s');
Kp = diag(kp);
Ki = diag(ki);
Kd = diag(kd);

K = Kp + Ki/s + Kd*s/(tauroll*s + 1);











