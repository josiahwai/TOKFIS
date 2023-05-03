% CONTROL DESIGN PROCEDURE TUTORIAL - NSTXU
clear all; clc; close all


%% STEP 1: INITIAL SCOPING

% preliminaries: tokamak geometry and power supply info
% tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
tok_fn = 'nstxu_obj_config2016_6565.mat';
utok = load(tok_fn).tok_data_struct;
circ = nstxu_circ();              % defines power supply limits and circuit connections
tok  = connect_tok(utok, circ);   % connect coils to circuits



%% STEP 1A: BUILD THE VACUUM MODEL
vacsys = response_models([], tok, 0, 0, 'vacuum');



%% STEP 1B: DESIGN A COIL CURRENT CONTROLLER
%
% Description: design using one of these different strategies 
% (pid, lqr, lqi) and then test it using cc_analysis.m

Kpid = nstxu_coilcontrol_pid(tok);  
Klqr = nstxu_coilcontrol_lqr(tok, vacsys);
Klqi = nstxu_coilcontrol_lqi(tok, vacsys);

K = Kpid;
stepsize = 1e3;
tmax = 0.2;
coils2plot = {'OH', 'PF5'};  
cc_analysis(K, vacsys, tok, stepsize, tmax, coils2plot, circ) 


%% STEP 1C: EVALUATE FLUX AND FIELD TEMPORAL RESPONSES
%
%  (can use to evaluate response speed for VS, shielding times, etc)
%  (details on whether to lock other coils, use coil control or voltage...)

r = [1 1.2];
z = [0 0];
Cbr = gridresponse2pt(tok.rg, tok.zg, [tok.gbr2c tok.gbr2v], r, z);
Cbz = gridresponse2pt(tok.rg, tok.zg, [tok.gbz2c tok.gbz2v], r, z);
Cpsi = gridresponse2pt(tok.rg, tok.zg, [tok.mpc tok.mpv], r, z);
Cvessel = [zeros(tok.nv,tok.nc) eye(tok.nv)];


stepsize = 1e7;  % T/A --> Gs/kA
t = linspace(0, 0.03, 100);
plotit = 1;
coils = tok.ccnames; 

t = linspace(0, 0.005, 100);
[~, figh] = vac_step_response(vacsys, K, Cbr, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Br [Gs/kA]', 'fontsize', 16)


[~, figh] = vac_step_response(vacsys, K, Cbz, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Bz [Gs/kA]', 'fontsize', 16)


stepsize = 1e6;
t = linspace(0, 0.3, 200);
[~, figh] = vac_step_response(vacsys, K, Cpsi, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Flux [mWb/kA]', 'fontsize', 16)


stepsize = 1e3;
t = linspace(0, 0.3, 200);
[~, figh] = vac_step_response(vacsys, K, Cvessel, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Vessel currents [A]', 'fontsize', 16)


T = Cbr*ss(vacsys.amat, vacsys.bmat, eye(size(vacsys.amat)), 0);
[y,t] = step(T(:,2), 0.005);
plot(t,y)

%% STEP 1D: EVALUATE FLUX AND FIELD SPATIAL RESPONSES

vec = [0 0 0 1 0 -1 0 0]';  % PF3U - PF3L
plot_response(tok, vec)

vec = [0 0.3 0.4 0.9 0 -0.9 -0.4 -0.3];
plot_response(tok, vec)


vec = [0 0 0 0 1 0 0 0];
plot_response(tok, vec)


vec = [1 0 0 0 0 0 0 0]';
plot_response(tok, vec);


%% STEP 1E: ESTIMATE OPTIMAL CONTROL DIRECTIONS

% Br optimization for vertical
r = linspace(0.6, 1.4, 10);
z = ones(size(r)) * 0.8;
r = [r r]';
z = [z -z]';
targ     = struct;
targ.r   = r;
targ.z   = z;
targ.br  = ones(size(r)) * 1e-7;
wt       = struct;
wt.br    = ones(size(r)) * 1;
vec = optimize_vec(tok, targ, wt);
plot_response(tok, normc(vec))


% Bz optimization for radial
z = linspace(-1,1,20)';
r = ones(size(z)) * 1.1;
targ     = struct;
targ.r   = r;
targ.z   = z;
targ.bz  = ones(size(r)) * 1e-7;
wt       = struct;
wt.bz    = ones(size(r)) * 1;
vec = optimize_vec(tok, targ, wt);
plot_response(tok, normc(vec))


% Flux optimization for Ip
r = linspace(0.6, 1.2, 10);
z = linspace(-1, 1, 10);
[r,z] = meshgrid(r,z);
r = r(:);
z = z(:);
targ      = struct;
targ.r    = r;
targ.z    = z;
targ.psi  = ones(size(r)) * 1e-6;
wt       = struct;
wt.psi   = ones(size(r)) * 1;
wt.ic    = [1 1 1 1 10 1 1 1]';
vec = optimize_vec(tok, targ, wt);
plot_response(tok, normc(vec))


figure
bar(categorical(tok.ccnames), normc(vec))


%% STEP 2A: DESIGN EQUILIBRIA

% Here we are just loading from a pre-saved file
% to create different equilibria, see generate_eqs.m
eqs_fn = './eq/eqs6565.mat';   
eqs = load(eqs_fn).eqs;



%% STEP 2B: BUILD STATE-SPACE MODELS

iplcirc = 1; 
Rp = 5e-7;
plasma_model = {'gspert', 'rzrig', 'rig'};

models = {};
for i = 1:length(eqs)
  models(:,end+1) = response_models(eqs{i}, tok, iplcirc, Rp, plasma_model);
end
models = models(:);



%% STEP 3A: ANALYZE VERTICAL INSTABILITY

gamma_max = 300;
gamma_min = -inf;
[iuse, gamma] = sys_eig_check(models, gamma_max, gamma_min);
models = models(iuse);
nmodels = length(models);

i = 3*4+1;
visualize_eig_vec(models{i}.amat, models{i}.eq, tok, 1)
sgtitle(models{i}.plasma_model)

i = i+1;
visualize_eig_vec(models{i}.amat, models{i}.eq, tok, 1)
sgtitle(models{i}.plasma_model)


%% STEP 3B: ESTIMATE CONTROL GAINS

eq = eqs{4};
gamma = 70;   % [Hz]
dt = 1/gamma; % sec
vec = [0 0 0 1 0 -1 0 0]';  % vector direction for control (+PF3U, -PF3L)
[kp0, kd0] = estimate_vertical_gains(eq, tok, vec, dt);


%% STEP 3C: OPTIMIZE CONTROL GAINS AGAINST MULTIPLE MODELS

% perform a gridscan to find stable regions
delay = 1e-4;          % [sec] time delay for vertical control system
opt.mag_up = 2;        % upper grid scan limit (order of magnitude)
opt.mag_lo = -2;       % lower grid scan limit (order of magnitude)
opt.model_order = 20;  % how much to compress model 
figure
scatter(log10(kp0), log10(kd0), 100, 'r', 'filled')   % the initial controller estimate
for i = 1:nmodels
  disp(i)
  vs_gridscan(models{i}, vec, kp0, kd0, delay, opt);  % find stable regions
end

% This script is a demo on how to design a controller using systune (which
% uses H-infinity optimization under the hood)
iuse = [4 14 15 25];
[Kz, K] = vs_optimization(models(iuse));


%% STEP 3D: ADD VERTICAL CURRENT COMPENSATION
%
% On NSTXU the vertical stability coil set (primarily PF3) are full-sized
% poloidal field coil, so current compensation is not needed. This is 
% important for KSTAR/SPARC where the vertical coils are fast but with
% smaller current capacity. See KSTAR/SPARC tutorials for analysis.


%% STEP 3E: ANALYZE AND TUNE LINEAR SYSTEM PERFORMANCE

iuse = [1:9 13:27];
delay = 1e-4;
opt = struct;
opt.reduce_model = 0;
opt.model_order  = 20;   
opt.stepsize = 0.01;
opt.t = linspace(0, 0.3, 500);

% analyze the controller made by systune in previous step
vs_analysis(models(iuse), Kz, delay, tok, opt);


% analyze this controller
s = tf('s');
vec = [0 0 0 1 0 -1 0 0]';
Kz1 = vec * (kp0 + kd0*s);
vs_analysis(models(iuse), Kz1, delay, tok, opt);


% analyze this controller
vec = [0 0.4 0.4 0.8 0 -0.8 -0.4 -0.4]';
Kz2 = vec * (kp0 + kd0*s);
vs_analysis(models(iuse), Kz2, delay, tok, opt);

%% STEP 4A: DESIGN RADIAL AND IP CONTROLLERS

sys = models{4};

% estimate radial gain
T = 0.05;
vec = [0 0 0 0 1 0 0 0]';
Rp = 5e-7;
iplcirc = 1;
vacsys = response_models(eq, tok, iplcirc, Rp, 'vacuum');
P = ss(vacsys.amat, vacsys.bmat, sys.drcurdx, 0);
y = step(P*vec, T);
kp = 5 ./ y(end) / eq.cpasma;
disp(kp)

% estimate Ip gain
vec = [1 0 0 0 0 0 0 0]';
Cip = [zeros(1,tok.nc+tok.nv) 1];
P = ss(vacsys.amat, vacsys.bmat, Cip, 0);
y = step(P*vec, T);
kp = 5 ./ y(end);
disp(kp)


%% STEP 4B/4C: ANALYZE LINEAR SYSTEM PERFORMANCE

iuse = [2:7 14:20];

opt.rstep = 0.01;  % [m]
opt.zstep = 0.01;  % [m]
opt.ipstep = 1e4;  % [A]
opt.t = linspace(0, 0.3, 500);
delay = 1e-4;

s = tf('s');

kp = 8e-3;
ki = 1e-2;
kd = 1e-4;
vec = [0 0 0 0 1 0 0 0]';
Kr = vec * (kp + ki/s + kd*s);

kp = -0.1;
ki = -0.4;
kd = 0;
vec = [1 0 0 0 0 0 0 0]';
Kip = vec * (kp + ki/s + kd*s); 

kp = 0.0021;
ki = 0;
kd = 1.2e-4;
vec = [0 0.4 0.4 0.8 0 -0.8 -0.4 -0.4]';
Kz = vec * (kp + ki/s + kd*s); 

rzip_analysis(models(iuse), Kr, Kz, Kip, delay, tok, opt)


%% STEP 4C: ANALYZE LINEAR SYSTEM PERFORMANCE WITH SIMULINK
sys = models{2};
linearsim_rzip(sys, tok, circ)


%% STEP 5: CHOOSE SHAPE CONTROL DESIGN
% This tutorial will cover both output-feedback and MPC


%% STEP 6: SHAPE INVERSE CONTROL
shape_cinv_controlsim        % this is a control simulation based on the 




%% STEP 10: NONLINEAR SIMULATION









































































