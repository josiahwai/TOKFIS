clear all; clc; close all

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
[r, z, s] = shape_create(s, 200);


targets.rsep = r;
targets.zsep = z;
weights.sep = ones(size(targets.rsep)) * 10;

targets.cpasma = 8e5;
weights.cpasma = 1e-3;

targets.psibry = 0.3;
weights.psibry = 10;

targets.ic = zeros(tok.nc,1);
weights.ic = ones(size(targets.ic)) * 1e-5;

targets.li = 1.2;
weights.li = 10;

targets.betap = 0.4;
weights.betap = 10;

targets.rx = nan;
targets.zx = nan;
weights.x  = 10;





% function [spec, init, config, eq] = gsdesign_from_scratch(targets, weights, tok, circ)
% config = tok;
% config.constraints = 1;
% config.psikn = (0:config.nr-1)/(config.nr-1);
% spec.limits.ic = [circ.icmax(:) circ.icmin(:)];









