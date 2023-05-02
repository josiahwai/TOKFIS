
function u = zcontroller(r, rdot, y, ydot, tnow, CONFIG)

e = r - y;
ed = rdot - ydot;
ip = CONFIG.ip;

vec = [0 0.2 0.2 0.8 0 -0.8 -0.2 -0.2]';
kp = 0.0021;
kd = 1.2e-4;

u = vec * ip * (kp*e + kd*ed);




