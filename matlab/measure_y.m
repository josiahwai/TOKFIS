
function y = measure_y(data, iy, targ, tok)


% initialize defaults
defaultvars = {'rxlo', 'zxlo', 'rxup', 'zxup', 'rxtarglo', ...
  'zxtarglo', 'rxtargup', 'zxtargup', 'rcp', 'zcp', 'rbdef', 'zbdef', ...
  'rtouch', 'ztouch', 'rtargtouch', 'ztargtouch'};
default = struct;
for x = defaultvars
  default.(x{:}) = nan;
end
targ = copyfields(default, targ, [], 1);

% 
psizr = reshape(data(iy.psizr), tok.nz, tok.nr);


% currents
y.ip = data(iy.cpasma);
y.ic = data(iy.ic);

% flux at control points
rcp = min(max(targ.rcp,tok.rg(3)), tok.rg(end-2));
zcp = min(max(targ.zcp,tok.zg(3)), tok.zg(end-2));

y.psicp = bicubicHermite(tok.rg, tok.zg, psizr, rcp(:), zcp(:));
y.psibry = data(iy.psibry);

% flux at target touch/x-points
[y.psixtarglo, y.psixtarglo_r, y.psixtarglo_z] = bicubicHermite(tok.rg, tok.zg, psizr, targ.rxtarglo, targ.zxtarglo);
[y.psixtargup, y.psixtargup_r, y.psixtargup_z] = bicubicHermite(tok.rg, tok.zg, psizr, targ.rxtargup, targ.zxtargup);
[y.psixlo, y.psixlo_r, y.psixlo_z] = bicubicHermite(tok.rg, tok.zg, psizr, targ.rxlo, targ.zxlo);
[y.psixup, y.psixup_r, y.psixup_z] = bicubicHermite(tok.rg, tok.zg, psizr, targ.rxup, targ.zxup);
y.psitouch = bicubicHermite(tok.rg, tok.zg, psizr, targ.rtouch, targ.ztouch);


% control points vs target touch/x-points
ONE = ones(length(rcp), 1);
y.diff_psicp_psixtarglo = y.psicp - ONE * y.psixtarglo;
y.diff_psicp_psixtargup = y.psicp - ONE * y.psixtargup;
y.diff_psicp_psibry     = y.psicp - ONE * y.psibry;
y.diff_psicp_psitouch   = y.psicp - ONE * y.psitouch;
y.diff_psicp_psixlo     = y.psicp - ONE * y.psixlo;




























