clc; close all

% Choose a model and equilibrium to use for this shape control simulation
i = 4; 
sys = models{i};
eq = sys.eq;

%% STEP 6A: DEFINE SHAPE VARIABLES TO CONTROL

% define isoflux control segments
opt = struct;
opt.rc = 0.85;
opt.zc = 0;
opt.a = 0.3; 
opt.b = 0.4;
opt.plotit = 0;
nsegs = 20;
segs = gensegs(nsegs, opt);


% Define target shape. In this example, we use the old shape and modify
% the elongation and triangularity
s = shape_analysis(eq.rbbbs, eq.zbbbs);
s.elong = 1.7;                              % modify elongation
s.triu = 0.4;                               % modify upper triangularity
[r, z] = shape_edit(eq.rbbbs, eq.zbbbs, s);
[zx,i] = min(z);
rx = r(i);
[r, z] = seg_intersections(r, z, segs, 0);

 
% Define which shape variables to control. (These are the only
% values that are explicitly controlled)
rzip_fds = {'rcur', 'zcur', 'cpasma'};
shape_fds = {'diff_psicp_psixlo', 'psixtarglo_r', 'psixtarglo_z'};
y_fds = [rzip_fds shape_fds];


% Specify target values
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
targ.rbdef = eq.rbdef;
targ.zbdef = eq.zbdef;

ref = struct2vec(targ, y_fds);    % shape control target in vector form


%% Define the coil current, R, Z, and Ip controllers

s = tf('s');

% Coil current controller. This will be used as part of the shape controller
kp = [4.3886    0.2447    0.2373    0.6230    1.4965    0.6230    0.2373    0.2447]';
ki = [87.7717    4.8936    4.7461   12.4610   29.9297   12.4610    4.7461    4.8936]';
kd = [0 0 0 0 0.1 0 0 0];
Kcoil = diag(kp) + diag(ki)/s + diag(kd)*s;

% Radial controller
kp = 8e-3;
ki = 1e-2;
kd = 1e-4;
vec = [0 0 0 0 1 0 0 0]';
Kr = vec * ip * (kp + ki/s + kd*s);


% Vertical controller
kp = 0.0021;
ki = 0;
kd = 1.2e-4;
vec = [0 0.2 0.2 0.8 0 -0.8 -0.2 -0.2]';
Kz = vec * ip * (kp + ki/s + kd*s); 


% Plasma current controller
kp = -0.1;
ki = -0.4;
kd = 0;
vec = [1 0 0 0 0 0 0 0]';
Kip = vec * (kp + ki/s + kd*s); 


%% DESIGN SHAPE CONTROLLER (STEP 6B and 6C)
% Includes step 6B: analyze and apply vertical decoupling
% Includes step 6C: compute weighted inverse plasma response


% Shape controller
dpsizrdx = tok.mpc;                        % use vacuum model for coil flux response
cdata = build_cdata(tok.mpc, tok, targ);   % build the C matrix 
C = struct2vec(cdata, shape_fds);
 
% define weights on the outputs
wy.psixlo            = 0;                              % flux at x-point
wy.diff_psicp_psixlo = ones(length(targ.rcp),1) * 5e9; % flux at control points minus flux at x-points
wy.psixtarglo_r      = 1e9;                            % flux gradient at x-point
wy.psixtarglo_z      = 1e9;                            % flux gradient at x-point

Wy = diag(struct2vec(wy,shape_fds));   % matrix form of weights

Wx = diag([1e3 1 1 1 1 1 1 1]);        % regularization weight on coil currents

% this block designs the regularization penalty to have low penalty on 
% (PF3U+PF3L) direction but high penalty on (PF3U-PF3L) direction
if 1                       
  vec1 = [0 0 0 1 0 +1 0 0];  % PF3U + PF3L
  vec2 = [0 0 0 1 0 -1 0 0];  % PF3U - PF3L  
  Dx = [eye(8); vec1; vec2];
  wx = [1e3 1 1 0 1 0 1 1 1 1e2];
  Wx = Dx' * diag(wx) * Dx;
end

H = C'*Wy*C + Wx;
f = -C'*Wy;
Cinv = -H\f;           % the weighted inverse plasma response

Kshape = Kcoil * Cinv; % shape controller is combo of coil current + shape mapping

%% Define simulation initial conditions and parameters

% power supply delay 
delay = 1e-4;  % sec
psdelay = pade(exp(-delay*s), 4);


% compute the full plasma flux response
mpp = unwrap_mpp(tok.mpp, tok.nz, tok.nr);   
dpsizrdx = mpp*sys.dpcurrtdx + [tok.mpc tok.mpv zeros(tok.nz*tok.nr,1)];   

response = build_cdata(dpsizrdx, tok, targ);
response.rcur = sys.drcurdx;
response.zcur = sys.dzcurdx;

% notation: 
%   x = coil currents + vessel currents + plasma current
%   y = outputs that are explicitly used for feedback (defined by y_fds)
%   z = all outputs that are measured

% define dydx for simulation
dydx = struct2vec(response, y_fds);
iy = get_struct_indices(response, y_fds);   % defines indices related to y

% define dzdx
z_fds = [y_fds tok.ccnames' {'psizr', 'psibry'} ];
dzdx = struct2vec(response, z_fds);
iz = get_struct_indices(response, z_fds);   % defines indices related to z


% define initial conditions

% x0: initial state
x0 = [eq.ic; eq.iv; eq.cpasma];             

% y0: initial feedback outputs
psicp = bicubicHermite(tok.rg, tok.zg, eq.psizr, targ.rcp, targ.zcp);
[~, eq.psixtarglo_r, eq.psixtarglo_z] = bicubicHermite(tok.rg, tok.zg, eq.psizr, targ.rxtarglo, targ.zxtarglo);

y0data = eq;
y0data.psizr = eq.psizr(:);
y0data.diff_psicp_psixlo = psicp - eq.psibry;   
y0 = struct2vec(y0data, y_fds);             

% z0: initial measured outputs
z0data = y0data;
for i = 1:tok.nc
  coil = tok.ccnames{i};
  z0data.(coil) = eq.ic(i);
end
z0 = struct2vec(z0data, z_fds);

%% STEP 6D: ANALYZE AND TUNE LINEAR SYSTEM PERFORMANCE 

% combine all controllers
K = [Kr Kz Kip Kshape];            

% plant model
P = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);  

% closed loop system
CL = feedback(P*psdelay*K, dydx);
[~,CL] = isproper(CL);

% simulate the system
n = 300;
t = linspace(0, 0.3, n);
[dx, t] = step(CL*(ref-y0), t);


% Reconstruct and plot outputs
% reconstruct the outputs
z = struct;
z.Time = t;
z.Data = z0 + dzdx * dx';

% plot output signals
ts = y2ts(z, iz);
fig = plot_structts(ts, tok.ccnames, 3);  % plot coil currents
fig.Position = [505 191 1069 709];
fig = plot_structts(ts, y_fds, 3);        % plot controlled outputs
fig.Position = [505 191 1069 709];

% plot flux 
fig = summary_soln_plot(z, iz, tok, targ);


% frequency analysis of the system
L = dydx*P*psdelay*K;        % loop gain r-->y

% nyquist plots
labels = {'rcur', 'zcur', 'cpasma', 'diff_psicp_psixlo1'};
idx = [iy.rcur iy.zcur iy.cpasma iy.diff_psicp_psixlo(1)];
h = plot_nyquists(L, idx, labels);








































