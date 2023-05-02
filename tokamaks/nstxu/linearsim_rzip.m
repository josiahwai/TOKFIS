% REQUIRED: 
%   r 
%   sys
%   PS.delay
%   PS.vmax
%   PS.vmin
%   dydx
%   x0
%   y0
%   CONFIG.controlfun_unicode
%   CONFIG.nu
%   simset.max_step_size
%   simset.stop_time

function out = linearsim_rzip(sys, tok, circ)

[nx,nu] = size(sys.bmat);
ny = 3;

PS.vmax = circ.vmax * inf;
PS.vmin = circ.vmin * inf;
PS.delay = 1e-4;

x0 = zeros(nx,1);
y0 = zeros(ny,1);

dydx = [sys.drcurdx; sys.dzcurdx; [zeros(1,tok.nc+tok.nv) 1]];

simset.max_step_size = 1e-3;
simset.stop_time = 0.3;

n = 100;
t = linspace(0, 0.3, n)';
rdata = zeros(n, 3);
rdata(t>0.001, 1) = 0.01;
r = timeseries(rdata, t);

controlfun_name = 'rzip_controller';

CONFIG = struct;
CONFIG.controlfun_unicode = double(controlfun_name); % pass the name of the controlfun, encoded as unicode vector
CONFIG.ip = sys.eq.cpasma;
CONFIG.iy = struct('rcur', 1, 'zcur', 2, 'cpasma', 3);
CONFIG.nu = 8;


% run model
out = run_slx_model('linearsim_slx', sys, PS, x0, y0, dydx, CONFIG, simset);
out.y.Data = squeeze(out.y.Data)';


% plot outputs
y = y2ts(out.y, CONFIG.iy);

fig = figure;
fig.Position = [680 296 485 681];
subplot(311)
plot(y.rcur)
subplot(312)
plot(y.zcur)
subplot(313)
plot(y.cpasma)























































