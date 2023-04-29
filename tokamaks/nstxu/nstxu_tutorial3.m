% CONTROL DESIGN PROCEDURE TUTORIAL - NSTXU
clear all; clc; close all


%% STEP 1: INITIAL SCOPING

% preliminaries: tokamak geometry and power supply info
tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
utok = load(tok_fn).tok_data_struct;
circ = nstxu_circ();              % defines power supply limits and circuit connections
tok  = connect_tok(utok, circ);   % connect coils to circuits



%% STEP 1A: BUILD THE VACUUM MODEL
sys = response_models([], tok, 0, 0, 'vacuum').vacuum;



%% STEP 1B: DESIGN A COIL CURRENT CONTROLLER
%
% Description: design using one of these different strategies 
% (pid, lqr, lqi) and then test it using cc_analysis.m

Kpid = nstxu_coilcontrol_pid(tok);  
Klqr = nstxu_coilcontrol_lqr(tok, sys);
Klqi = nstxu_coilcontrol_lqi(tok, sys);

K = Kpid;
stepsize = 1e3;
tmax = 0.2;
coils2plot = {'OH', 'PF5'};  
cc_analysis(K, sys, tok, stepsize, tmax, coils2plot, circ) 


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
[~, figh] = vac_step_response(sys, K, Cbr, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Br [Gs/kA]', 'fontsize', 16)


[~, figh] = vac_step_response(sys, K, Cbz, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Bz [Gs/kA]', 'fontsize', 16)


stepsize = 1e6;
t = linspace(0, 0.3, 200);
[~, figh] = vac_step_response(sys, K, Cpsi, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Flux [mWb/kA]', 'fontsize', 16)


stepsize = 1e3;
t = linspace(0, 0.3, 200);
[~, figh] = vac_step_response(sys, K, Cvessel, tok, circ, stepsize, t, coils, plotit);
tabfun(figh, @ylabel, 'Vessel currents [A]', 'fontsize', 16)


T = Cbr*ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);
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

figure
bar(categorical(tok.ccnames), normc(vec))

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

sys = {};
for i = 1:length(eqs)
  sys(:,end+1) = response_models(eqs{i}, tok, iplcirc, Rp, plasma_model);
end



%% STEP 3A: ANALYZE VERTICAL INSTABILITY

sys = sys(:);
[iuse, gamma] = sys_eig_check(sys, 300);
sys = sys(iuse);

i = 3*4+1;
visualize_eig_vec(sys{i}.amat, sys{i}.eq, tok, 1)
sgtitle(sys{i}.plasma_model)

i = i+1;
visualize_eig_vec(sys{i}.amat, sys{i}.eq, tok, 1)
sgtitle(sys{i}.plasma_model)


%% STEP 3B: ESTIMATE CONTROL GAINS

dz = 0.03;  % m
gamma = 70; % Hz

T = 1 / gamma;

br = tok.gbr2c*eq.ic + tok.gbr2v*eq.iv;
br = reshape(br, tok.nz, tok.nr);
dr = mean(diff(tok.rg));
dz = mean(diff(tok.zg));
[~, dbrdz] = gradient(br, dr, dz);

r = linspace(0.6, 1.2, 10);
z = linspace(-1, 1, 10);
[r,z] = meshgrid(r,z);

i = tok.zgg < 1 & tok.zgg > -1 & tok.rgg > 0.6 & tok.rgg < 1.2;
dbrdz = mean(dbrdz(i));

sys = response_models([], tok, 0, 0, 'vacuum');
sys = sys{1};

r = 1;
z = 0;
Cbr = gridresponse2pt(tok.rg, tok.zg, [tok.gbr2c tok.gbr2v], r, z);

OL = Cbr*ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);
[y,t] = step(OL(:,4), T);

plot(t,y)
dbrdv = y(end);

kp0 = 0.1 * dbrdz / dbrdv / eq.cpasma;
kd0 = 0.5 * dbrdz * T / dbrdv / eq.cpasma;


%%
close all

figure
for i = 1:length(syss)
  i
  sys = syss{i};
  kp0 = 0.002;
  kd0 = 1.6e-4;
  vec = [0 0 0 1 0 -1 0 0]';

  [gamma, kps, kds] = vs_gridscan(sys, vec, kp0, kd0);
end

%%
% REPEAT WITH TIME DELAY
i = 1;
sys = syss{i};
tau = 2e-3;
[gamma, kps, kds] = vs_gridscan_timedelay(sys, vec, kp0, kd0, tau);

scatter(log10(kp0), log10(kd0), 100,  'k', 'filled')

%%
% close all

figure
for i = 1:length(syss)
  i
  sys = syss{i};
  kp0 = 0.002;
  kd0 = 1.6e-4;
  vec = [0 0 0 1 0 -1 0 0]';
  tau = 1e-3;
  opt.model_order = 20;
  
  [gamma, kps, kds] = vs_gridscan_timedelay(sys, vec, kp0, kd0, tau, opt);
end

scatter(log10(kp0), log10(kd0), 100,  'k', 'filled')


%%
close all
kp = 0.008;
kd = .0002;

s = tf('s');
P = ss(sys.amat, sys.bmat*vec, sys.dzcurdx, 0);
K = ip*kp + ip*kd*s;
K = K * exp(-0.003*s);

CL = feedback(P*K, 1);
max(real(pole(pade(CL,2))))


hold on
step(CL, 0.1)


  
%%
% 
% s = struct('rsurf', 0.9, 'aminor', 0.55, 'elong', 1.5, 'c_xplo', 0);
% [r,z] = shape_create(s, 100);
% clf; plot_lim(tok); plot(r,z)



uigrid_plot(tok, x, v)

% eqs_fn = './eq/eqs6565.mat';
% 
% 
% 
% 
% 
% % LOAD SOME EQUILIBRIA AND GEOMETRY
% % (to create different equilibria, see generate_eqs.m)
% tok = load(tok_fn).tok_data_struct;
% eqs = load(eqs_fn).eqs;
% eq = eqs{3};
% 
% % CONNECT COILS-->CIRCUITS
% circ = nstxu_circ();
% tok = connect_tok(tok, circ);          
% 
% 
% % BUILD THE VACUUM MODEL
% sys = response_models([], tok, 0, 0, 'vacuum');
% sys = sys.vacuum;



% DESIGN A COIL CURRENT CONTROLLER
%
% K = nstxu_coilcontrol_pid();
% K = nstxu_coilcontrol_lqr();
% K = nstxu_coilcontrol_lqi();
% vmax = 
% vmin = 
% coil_current_analysis(sys, K, 'CS1U', vmax, vmin)



% EVALUATE FLUX AND FIELD TEMPORAL RESPONSES
%  (can use to evaluate response speed for VS, shielding times, etc)
%  (details on whether to lock other coils, use coil control or voltage...)
%
% r = 1.2;
% z = 0;
% vec = 
% appvolt = 1e3;
% tok_volt_step_response(sys, 'Br', r, z, vec, appvolt)
%
% tok_ref_current_step_response(sys, 'Br', r, z, vec, )



% EVALUATE FLUX AND FIELD SPATIAL RESPONSE
%
% dI = 
% dI(iy.OH) = 
% visualize_response(sys, 'Br', dI)



% OPTIMIZE FOR FLUX AND FIELD SPATIAL RESPONSE OVER A PARTICULAR AREA
% 
% targets = 
% weights = 
% dx = optimize_dx(sys, targets, weights)



% BUILD PLASMA RESPONSE MODELS
%
% iplcirc = 1; 
% Rp = 5e-7;
% models = response_models(eq, tok, iplcirc, Rp, plasma_model);



% VISUALIZE/FAMILIARIZE A BIT WITH RESPONSE MODEL
%
% visualize_eig_vec(sys, 1)
% visualize_response(sys, 'Br', dI)
% tok_ref_current_step_response(sys, 'Br', r, z, vec, )



% VERTICAL CONTROL DEVELOPMENT
%
% iplcirc = 0;
% models = response_models(eq, tok, iplcirc, Rp, plasma_model);
% z = 0;
% r = 1.2;
% tok_ref_current_step_response(sys, 'Br', r, z, vec, )
%
% [kp0, kd0] = zn_tuning(...)
%
% [kp0, kd0] = estimate_zcontrol_params(sys)
% [kp, kd] = gridsearch_stable_zcontrol(sys, kp0, kd0)
% tauroll = 1e-6; 
% Kz = kp + kd*s/(tauroll*s + 1)
% fluxloop_vec = dz2dfluxloop(sys)
% zcontrol_nyquist(sys, Kz)
% zcontrol_bode(sys, Kz)
% zcontrol_step_response(sys, Kz)
% zcontrol_disturb_response(sys, Kz)
%
% now apply to multiple models
% sys = {sys1, sys2, sys3...}
% zcontrol_step_response(sys, Kz)
% zcontrol_montecarlo_hinfstruct(sys, Kz0)
% [kp, ki, kd] = control_fminsearch(sys1, sys2, ...)


 
% RADIAL CONTROL DEVELOPMENT
% 
% vec = 
% [kp, ki, kd] = zn_tuning(sys ...)
% [kp, ki, kd] = robusttune(sys1, sys2, ...)
% [kp, ki, kd] = control_fminsearch(sys1, sys2, ...)
% linear_sim(?)



% IP CONTROL DEVELOPMENT
% 
% vec = 
% [kp, ki, kd] = zn_tuning(sys ...)
% [kp, ki, kd] = robusttune(sys1, sys2, ...)
% [kp, ki, kd] = control_fminsearch(sys1, sys2, ...)
% linear sim(?)



% CINV-BASED SHAPE CONTROL DEVELOPMENT
%
% cdata = build_cdata(dpsizrdx)
% Kshape = cmat_inversion(cdata, wts, targets, fds2control ...)
% control_freq_analysis(sys, Kshape, Kcoil, Kz, Kr, Kip ...)
% linear_sim(sys, K...)
%
% TRY AGAIN WITH VERTICAL COMPENSATION:
% cdata = build_cdata(dpsizrdx)
% u = shapecontrol(y, c.iy, tnow, config)
% linearsim(@shapecontrol, sys, eq, K ...)
%
% FIND RHP ZERO DIRECTION AND PENALIZE IT
% rhp_zero_analysis(...)
% wt = 
% K = cmat_inversion(wt, ....)
% coil_current_step_response(vec1, ...
% coil_current_step_response(vec2, ...
%
% Can cross-test with other equilibria, sys (i.e.)
% linearsim(sys2, eq, K ...)



% MODEL PREDICTIVE CONTROL 
%
% define_lotsa_stuff(...)
% control_freq_analysis(...)
% linear_sim(...)
% also penalize RHP zero direction



% NONLINEAR SIMULATIONS - LOTS OF EXAMPLES



































































