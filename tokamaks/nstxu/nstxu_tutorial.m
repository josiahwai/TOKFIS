% WORKFLOW TUTORIAL FOR DESIGNING A SHAPE CONTROLLER - NSTXU
clear all; clc; close all

tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
eqs_fn = './eq/eqs6565.mat';

% LOAD SOME EQUILIBRIA AND GEOMETRY
% (to create different equilibria, see generate_eqs.m)
tok = load(tok_fn).tok_data_struct;
eqs = load(eqs_fn).eqs;
eq = eqs{3};

% CONNECT COILS-->CIRCUITS
circ = nstxu_circ();
tok = connect_tok(tok, circ);          


% BUILD THE VACUUM MODEL
sys = response_models([], tok, 0, 0, 'vacuum');
sys = sys.vacuum;



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



































































