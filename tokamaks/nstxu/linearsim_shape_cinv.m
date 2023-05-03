clear; clc; close all

load('x')
i = 4;
sys = models{i};
eq = sys.eq;

% circ = nstxu_circ();             
% tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
% tok = load(tok_fn).tok_data_struct;
% tok  = connect_tok(tok, circ);   
% 
% eqs_fn = './eq/eqs3333.mat';   
% eq = load(eqs_fn).eqs{2};
% 
% iplcirc = 1; 
% Rp = 5e-7;
% plasma_model = 'gspert';
% sys = response_models(eq, tok, iplcirc, Rp, plasma_model);

simset.max_step_size = 1e-3;
simset.stop_time = 0.3;
simset.decimation = 1;



%%

% define control segments
opt = struct;
opt.rc = 0.85;
opt.zc = 0;
opt.a = 0.3; 
opt.b = 0.4;
opt.plotit = 0;
nsegs = 20;
segs = gensegs(nsegs, opt);


% define target shape
s = shape_analysis(eq.rbbbs, eq.zbbbs);
s.elong = 1.7;
s.triu = 0.4;
[r, z] = shape_edit(eq.rbbbs, eq.zbbbs, s);
[zx,i] = min(z);
rx = r(i);
[r, z] = seg_intersections(r, z, segs, 0);

targ.rcp = r;
targ.zcp = z;
targ.rxlo = rx;
targ.zxlo = zx;
targ.rxtarglo = rx;
targ.zxtarglo = zx;
targ.rcur = eq.rcur;
targ.zcur = eq.zcur;
targ.cpasma  = eq.cpasma;
targ.diff_psicp_psixlo = zeros(nsegs,1);
targ.psixtarglo_r = 0;
targ.psixtarglo_z = 0;
targ.psixlo = 0;

mpp = unwrap_mpp(tok.mpp, tok.nz, tok.nr);   

eqx = targ;
eqx.rbdef = eq.rbdef;
eqx.zbdef = eq.zbdef;

dpsizrdx = mpp*sys.dpcurrtdx + [tok.mpc tok.mpv zeros(tok.nz*tok.nr,1)];

cdata = build_cdata(dpsizrdx, tok, eqx);
cdata.rcur = sys.drcurdx;
cdata.zcur = sys.dzcurdx;

yfds = {'rcur', 'zcur', 'cpasma', 'ic', 'psizr', 'psibry'};
iy = get_struct_indices(cdata, yfds);
dydx = struct2vec(cdata, yfds);

[nx,nu] = size(sys.bmat);
ny = size(dydx,1);

x0 = [eq.ic; eq.iv; eq.cpasma];

y0 = eq;
y0.psizr = y0.psizr(:);
y0 = struct2vec(y0, yfds);

PS.vmax = circ.vmax * inf;
PS.vmin = circ.vmin * inf;
PS.delay = 1e-4;


controlfun_name = 'shape_cinv_controller';

t = linspace(0, 0.3, 100)';
r = struct;
r.Time = t;
r.Data = ones(length(t),1) * y0';
% r.Data(t>0.1, 2) = 0.01;

% r.Data = r.Data(:,1:3);

CONFIG = struct;
CONFIG.controlfun_unicode = double(controlfun_name); % pass the name of the controlfun, encoded as unicode vector
CONFIG.ip = sys.eq.cpasma;
CONFIG.iy = iy;
CONFIG.nu = 8;
CONFIG.targ = targ;
CONFIG.dpsizrdx = tok.mpc;
CONFIG.tok = copyfields(struct, tok, {'rg','zg','nr','nz', 'nc', 'nv'}, 1);
CONFIG.r = r;
CONFIG = purge_struct(purge_struct(CONFIG));

%%
clc; close all

% simulate model
clear(controlfun_name)
% out = run_slx_model('linearsim_slx', sys, PS, x0, y0, dydx, CONFIG, simset);
% out.y.Data = squeeze(out.y.Data)';

tic
out = sim('linearsim_slx', 0.3);
toc


% read outputs
for x = tok.ccnames(:)'
  coil = x{:};
  iy.(coil) = iy.ic(circ.iy.(coil));
end
ts = y2ts(out.y, iy);


% plot outputs
fds2plot = [{'rcur', 'zcur', 'cpasma'} tok.ccnames'];
fig = plot_structts(ts, fds2plot, 3);
fig.Position = [505 191 1069 709];


fig = summary_soln_plot(out.y, iy, tok, targ);



















































