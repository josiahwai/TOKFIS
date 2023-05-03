clear; clc; close all

load('x')
i = 4;
sys = models{i};
eq = sys.eq;

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

targ.rbdef = eq.rbdef;
targ.zbdef = eq.zbdef;

dpsizrdx = mpp*sys.dpcurrtdx + [tok.mpc tok.mpv zeros(tok.nz*tok.nr,1)];

cdata = build_cdata(dpsizrdx, tok, targ);
cdata.rcur = sys.drcurdx;
cdata.zcur = sys.dzcurdx;

rzip_fds = {'rcur', 'zcur', 'cpasma'};
shape_fds = {'diff_psicp_psixlo', 'psixtarglo_r', 'psixtarglo_z'};
y_fds = [rzip_fds shape_fds];

iy = get_struct_indices(cdata, y_fds);
dydx = struct2vec(cdata, y_fds);

[nx,nu] = size(sys.bmat);
ny = size(dydx,1);

x0 = [eq.ic; eq.iv; eq.cpasma];


psicp = bicubicHermite(tok.rg, tok.zg, eq.psizr, targ.rcp, targ.zcp);
[~, eq.psixtarglo_r, eq.psixtarglo_z] = bicubicHermite(tok.rg, tok.zg, eq.psizr, targ.rxtarglo, targ.zxtarglo);

y0 = eq;
y0.diff_psicp_psixlo = psicp - eq.psibry;
y0 = struct2vec(y0, y_fds);

ref = struct2vec(targ, y_fds);
%%

ip = eq.cpasma;

Kpid = nstxu_coilcontrol_pid(tok);  

s = tf('s');

kp = 8e-3;
ki = 1e-2;
kd = 1e-4;
vec = [0 0 0 0 1 0 0 0]';
Kr = vec * ip * (kp + ki/s + kd*s);

kp = -0.1;
ki = -0.4;
kd = 0;
vec = [1 0 0 0 0 0 0 0]';
Kip = vec * (kp + ki/s + kd*s); 

kp = 0.0021;
ki = 0;
kd = 1.2e-4;
vec = [0 0.4 0.4 0.8 0 -0.8 -0.4 -0.4]';
Kz = vec * ip * (kp + ki/s + kd*s); 

%%

cdata = build_cdata(tok.mpc, tok, targ);
C = struct2vec(cdata, shape_fds);
 
wy.psixlo            = 0;
wy.diff_psicp_psixlo = ones(length(targ.rcp),1) * 3e10;
wy.psixtarglo_r      = 6e9;
wy.psixtarglo_z      = 6e9;

Wy = diag(struct2vec(wy,shape_fds));
Wx = diag([1e3 1 1 1 1 1 1 1]);

H = C'*Wy*C + Wx;
f = -C'*Wy;
Cinv = -H\f;

Kshape = Kpid * Cinv;

K = [Kr Kz Kip Kshape];


P = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);

CL = feedback(P*K, dydx);
[~,CL] = isproper(CL);

%%

K = [Kr Kz Kip Kshape*0];
ref = y0;
ref(2) = y0(2) + 0.01;

% n = 300;
% t = linspace(0, 0.3, n);
% u = (ref - y0) * ones(1,n);
% dx = lsim(CL, u, t);

[dx, t] = step(CL*(ref-y0), 0.3);

y.Time = t;
y.Data = dydx*dx';


% plot outputs
ts = y2ts(y, iy);
fig = plot_structts(ts, y_fds, 3);
fig.Position = [505 191 1069 709];

%%

% % read outputs
% for x = tok.ccnames(:)'
%   coil = x{:};
%   iy.(coil) = iy.ic(circ.iy.(coil));
% end
% ts = y2ts(out.y, iy);
% 
% 
% % plot outputs
% fds2plot = [{'rcur', 'zcur', 'cpasma'} tok.ccnames'];
% fig = plot_structts(ts, fds2plot, 3);
% fig.Position = [505 191 1069 709];
% 
% 
% fig = summary_soln_plot(out.y, iy, tok, targ);



















































