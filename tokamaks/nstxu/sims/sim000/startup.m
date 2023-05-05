clear; clc; close all

setenv('GATOOLS_ROOT', '/Users/jwai/Research/toksys/')
gatools_root = getenv('GATOOLS_ROOT');
run([gatools_root 'startups/toksys_startup'])

%% load equilibrium and tok
tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
eqs_fn = 'eqs3333.mat';

circ = nstxu_circ();
eq = load(eqs_fn).eqs{2};
tok = load(tok_fn).tok_data_struct;
tok = connect_tok(tok, circ);

init = eq;
%% GSEVOLVE SETUP

s = struct;     % will hold all the simulation settings
s.vacuum     = 0;
s.start_time = 0;
s.stop_time  = 0.95;
s.step_time  = 1e-6;
s.dt2save    = 50e-6;
s.decimation = round(s.dt2save/s.step_time);
s.init       = init;
s.tok        = tok;

t = linspace(0,1);
s.betap = timeseries(init.betap, t);
s.li    = timeseries(init.li, t);
s.Rp    = timeseries(5e-7, t);

% clear gsupdate cache
clear gsupdate

Tokamak = nstxu_Tokamak(s);
Actuators = nstxu_Actuators();
c = Tokamak.c;
e = Tokamak.e;




















































