% Use LQI to design the coil current controller. Similar to LQR but
% also has control on the integral of errors. 

function [K, Kp, Ki] = nstxu_coilcontrol_lqi(tok, sys)

A = sys.amat;
B = sys.bmat;
C = eye(tok.nc, tok.nc+tok.nv);
P = ss(A,B,C,0);

L = diag(tok.mcc);   
wt.ic = (L*180).^2;            % weight on coil currents
wt.iv = zeros(tok.nv,1);       % weight on vessel currents is zero
wt.ic_integral = wt.ic * 1000;  % weight on integral of coil currents


% LQI weighting matrices
Q = diag([wt.ic; wt.iv; wt.ic_integral]);
R = eye(tok.nc);


% solve LQI
tmp = lqi(P, Q, R);
Kp = tmp(:,1:tok.nc);       % proportional gains, dropping the vessel terms 
Ki = -tmp(:,end-tok.nc+1:end);  % integral gains


% make controller
s = tf('s');
K = Kp + Ki/s;



