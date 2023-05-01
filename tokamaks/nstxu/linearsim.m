clc; close all

i = 4;

sys = models{i};
[nx,nu] = size(sys.bmat);

PS.vmax = circ.vmax * inf;
PS.vmin = circ.vmin * inf;
PS.delay = 1e-3;


x0 = zeros(nx,1);
y0 = 0;

dydx = sys.dzcurdx;
r = timeseries([0 0 1 1], [0 0.1 0.101 0.3]);


simset.max_step_size = 1e-3;
simset.stop_time = 0.3;


controlfun_name = 'zcontroller';




CONFIG = struct;

% pass the name of the controlfun, encoded as unicode vector
CONFIG.controlfun_unicode = double(controlfun_name);   


CONFIG.ip = sys.eq.cpasma;
CONFIG.nu = nu;

out = sim('linearsim_slx');

figure
plot(out.y)






















































