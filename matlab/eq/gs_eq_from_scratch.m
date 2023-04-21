clear all; clc; close all


function [spec, init, config, eq] = gsdesign_from_scratch(targ, tok, circ)


tok = load('nstxu_obj_2016_GSgrid33x33_npp4x4.mat').tok_data_struct;
circ = nstxu_circ();
tok = connect_tok(tok, circ);


s = struct;
s.rsurf = 0.93;
s.zsurf = -0.02;
s.aminor = 0.56;
s.elong = 1.8;
s.triu = 0.3;
s.tril = 0.6;
s.squo = -0.1;
s.sqlo = -0.2;
s.squi = -0.1;
s.sqli = -0.3;
s.c_xplo = 0.1;
s.c_xpup = 0;
[rsep, zsep, s] = shape_create(s, 200);


config = tok;
config.constraints = 1;
config.psikn = (0:config.nr-1)/(config.nr-1);


spec.targets.cpasma = 8e5;
spec.weights.cpasma = 1e-3;

spec.targets.psibry = 0.3;
spec.weights.psibry = 10;

spec.targets.ic = zeros(tok.nc,1);
spec.weights.ic = ones(size(spec.targets.ic)) * 1e-5;

spec.targets.li = 1.2;
spec.weights.li = 10;

spec.targets.betap = 0.4;
spec.weights.betap = 10;

spec.targets.rx = nan;
spec.targets.zx = nan;
spec.weights.x  = 10;

spec.limits.ic = [circ.icmax(:) circ.icmin(:)];









