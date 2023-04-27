% CONTROL DESIGN PROCEDURE TUTORIAL - NSTXU
clear all; clc; close all


%% STEP 1: INITIAL SCOPING

% preliminaries: tokamak geometry and power supply info
tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
tok = load(tok_fn).tok_data_struct;
circ = nstxu_circ();            % defines power supply limits and circuit connections
tok = connect_tok(tok, circ);   % connect coils to circuits



% STEP 1A: BUILD THE VACUUM MODEL
sys = response_models([], tok, 0, 0, 'vacuum').vacuum;



% STEP 1B: DESIGN A COIL CURRENT CONTROLLER
%
% Description: design using one of these different strategies 
% (pid, lqr, lqi) and then test it using cc_analysis.m

Kpid = nstxu_coilcontrol_pid(tok);  
Klqr = nstxu_coilcontrol_lqr(tok, sys);
Klqi = nstxu_coilcontrol_lqi(tok, sys);

% K = Klqr;
stepsize = 1e3;
tmax = 0.2;
% coils2plot = 'all';  
% cc_analysis(K, sys, tok, stepsize, tmax, coils2plot) 
% 
% K(2,1) = 0.2;
% K(8,1) = 0.2;
% coils2plot = 'OH';
% cc_analysis_vlim(K, sys, tok, stepsize, tmax, coils2plot, circ)  
% 
% K = Kpid;
% cc_analysis_vlim(K, sys, tok, stepsize, tmax, coils2plot, circ)  

[~,Kp,Ki] = nstxu_coilcontrol_lqi(tok, sys);
Kp(2,1) = 0.2;
Kp(8,1) = 0.2;
Ki(2,1) = 0;
Ki(8,1) = 0;
s = tf('s');
K = Kp + Ki/s;
coils2plot = {'OH', 'PF5'};
cc_analysis_vlim(K, sys, tok, stepsize, tmax, coils2plot, circ)  


%%

% STEP 1C: EVALUATE FLUX AND FIELD TEMPORAL RESPONSES
%  (can use to evaluate response speed for VS, shielding times, etc)
%  (details on whether to lock other coils, use coil control or voltage...)







[x,t] = coilsim()

A = sys.amat;
B = sys.bmat;
Cc = eye(tok.nc, tok.nc+tok.nv);
P = ss(A,B,Cc,0);              % plant model (gain from u to y)
G = ss(A,B,eye(size(A)), 0);  % plant model including vessels

K = ss(K);                % controller model
I = eye(size(P));         % identity matrix
L   = P*K;                % loop transfer function
Gry = feedback(P*K, I);   % closed-loop transfer function (gain from r to y)
Gre = feedback(I, P*K);   % sensitivity function (gain from r to e)
Gru = feedback(K, P);     % control effort (gain from r to u)


n = 400;
t = linspace(0, tmax, n);

simIn = Simulink.SimulationInput('coilsim');
simIn = setVariable(simIn,'G',G);
simIn = setVariable(simIn,'K',K);
simIn = setVariable(simIn,'circ',circ);
simIn = setVariable(simIn,'tmax',tmax);
simIn = setVariable(simIn, 'Cc', Cc);



  r = zeros(tok.nc,1);
  r(i) = stepsize;
  simIn = setVariable(simIn,'r',r);


% open loop plant including vessels 
% (transfer function from voltage to currents)
P = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);  


% closed loop plant including vessels
% (transfer function from current targ  current)
Cx = eye(tok.nc, tok.nc+tok.nv);
T = feedback(P*K, Cx);                               


r = 1;
z = 0;
Cbr = gridresponse2pt(tok.rg, tok.zg, [tok.gbr2c tok.gbr2v], r, z);
Cbz = gridresponse2pt(tok.rg, tok.zg, [tok.gbz2c tok.gbz2v], r, z);
Cpsi = gridresponse2pt(tok.rg, tok.zg, [tok.mpc tok.mpv], r, z);


[x,t] = step(T(:,4), 0.5);
plot(t,x)

figure
plot(t, x*Cbz')

[x,t] = step(P(:,4), 0.05);
figure
plot(t,x)

figure
plot(t, x*Cbr')



%%


% tok_volt_step_response(sys, 'Br', r, z, vec, appvolt)
% tok_ref_current_step_response(sys, 'Br', r, z, vec, )


%
% r = 1.2;
% z = 0;
% vec = 
% appvolt = 1e3;

r = 1;
z = 0;
gbr = [tok.gbr2c tok.gbr2v];
C = [];
for i = 1:size(gbr,2)
  x = reshape(gbr(:,i), tok.nz, tok.nr);
  C(i) = bicubicHermite(tok.rg, tok.zg, x, r, z);
end
Gbr = ss(sys.amat, sys.bmat, C, 0);
step(Gbr(8), 0.05)
grid on



%%


% 
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



































































