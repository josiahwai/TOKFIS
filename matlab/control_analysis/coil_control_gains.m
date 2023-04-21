% coil current control
% PID, LQR, LQI, dLQR 
% look at voltages and cross-disturbances for each for step responses
% vmax
% vmin

clear all; clc; close all
load('x')


%% Coil current control
sys = response_models(eq, tok, 0, 0, 'vacuum');
sys = sys.vacuum; 

% compress vessel elements
nvessmodes = 40;
sys = compress_sys_vessels(sys, nvessmodes, tok);

% build the plant state-space model
P = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);
[nx, nu] = size(P);



%% PID-based control

Kp_pid = diag([3 3 1 1 1 1 1 1 5 5 5 5 10 10 0.3 0.3 0.3 0.3 1]) * 4; 
Ki_pid = diag([1 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 0]);
Kd_pid = Kp_pid * 0;

tauroll = 1e-2;


%% LQR-based control

coil_wts = [2 2 1 1 1 1 1 1 1 1 1 1 2 2 0.1 0.1 0.1 0.1 0.1] * 200;
vessel_wts = zeros(1,nvessmodes);

Q = diag([coil_wts vessel_wts]);
R = eye(nu);

Kp_lqr = lqr(sys.amat, sys.bmat, Q, R);
Kp_lqr = Kp_lqr(:,1:tok.nc);


%% LQI-based control 
ny = tok.nc;

A = sys.amat;
B = sys.bmat;
C = eye(ny, nx);

coil_wts = [2 2 1 1 1 1 1 1 1 1 1 1 2 2 0.1 0.1 0.1 0.1 0.1] * 200;
vessel_wts = zeros(1, nvessmodes);
coil_integral_wts = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * 1e4;

Q = diag([coil_wts vessel_wts coil_integral_wts]);
R = eye(nu); 

tmp = ss(A,B,C,0);
Klqi = lqi(tmp, Q, R);

Kp_lqi = Klqi(:,1:tok.nc);
Ki_lqi = -Klqi(:,end-tok.nc+1:end);































