
function u = rzip_controller(r, rdot, y, ydot, tnow, CONFIG)

persistent ei tprev
if isempty(ei), ei = r*0; end
if isempty(tprev), tprev = tnow; end

e = r-y;                   % errors
ei = ei + e*(tnow-tprev);  % integral of errors
edot = rdot - ydot;

u = zeros(8,1);

ip = CONFIG.ip;
iy = CONFIG.iy;


% radial control
vec = [0 0 0 0 1 0 0 0]';
kp = 8e-3;
ki = 1e-2;
kd = 2e-4;
c = 0.5;
tmp = kp*e(iy.rcur) + ki*ei(iy.rcur) + kd * (c*rdot(iy.rcur) - ydot(iy.rcur));
u = u + vec * ip * tmp;

 
% vertical control
vec = [0 0.2 0.2 0.8 0 -0.8 -0.2 -0.2]';
kp = 0.0021;
ki = 0;
kd = 1.2e-4;
tmp = kp*e(iy.zcur) + ki*ei(iy.zcur) + kd*edot(iy.zcur);
u = u + vec * ip * tmp;


% Ip control
vec = [1 0 0 0 0 0 0 0]';
kp = -0.1;
ki = -0.4;
kd = -0.001;
c = 0.5;
tmp = kp*e(iy.cpasma) + ki*ei(iy.cpasma) + kd * (c*rdot(iy.cpasma) - ydot(iy.cpasma));
u = u + vec * tmp;


tprev = tnow;

