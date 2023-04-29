% Estimate vertical control gains for vertical proportional derivative controller
%
% Method: correlate dBr/dz and dBr/dt
%
%

function [kp, kd] = estimate_vertical_gains(eq, tok, vec, dt)


% calculate equilibrium dBr/dz on grid
br = tok.gbr2c*eq.ic + tok.gbr2v*eq.iv;
br = reshape(br, tok.nz, tok.nr);
dr = mean(diff(tok.rg));
dz = mean(diff(tok.zg));
[~, dbrdz] = gradient(br, dr, dz);


% compute an average dBr/dz value
rmin = quantile(tok.rg, 0.4);
rmax = quantile(tok.rg, 0.6);
zmin = quantile(tok.rg, 0.4);
zmax = quantile(tok.rg, 0.6);
i = tok.zgg < zmax & tok.zgg > zmin & tok.rgg > rmin & tok.rgg < rmax;

dbrdz = mean(abs(dbrdz(i))) * median(sign(dbrdz(tok.zgg>0)));



% compute a vacuum step response of radial field 
sys = response_models([], tok, 0, 0, 'vacuum');
sys = sys{1};

r = median(tok.rg);
z = median(tok.zg);
Cbr = gridresponse2pt(tok.rg, tok.zg, [tok.gbr2c tok.gbr2v], r, z);
P = ss(sys.amat, sys.bmat*vec, eye(size(sys.amat)), 0);
OL = Cbr*P;
y = step(OL, dt);

dbrdv = y(end);  % how much br has changed due to applied voltage after dt


% estimate the control gains

a = 0.2;   % assume 0.2 contribution from proportional, 0.8 from derivative
kp = a * dbrdz / dbrdv / eq.cpasma;
kd = (1-a) * dbrdz * dt / dbrdv / eq.cpasma;


























