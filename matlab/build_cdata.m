
function cdata = build_cdata(dpsidx, tok, eqx)

% initialize defaults
defaultvars = {'rxlo', 'zxlo', 'rxup', 'zxup', 'rxtarglo', ...
  'zxtarglo', 'rxtargup', 'zxtargup', 'rcp', 'zcp', 'rbdef', 'zbdef', ...
  'rtouch', 'ztouch', 'rtargtouch', 'ztargtouch'};
default = struct;
for x = defaultvars
  default.(x{:}) = nan;
end
eqx = copyfields(default, eqx, [], 1);

% response of conductors is identity
c.x = eye(tok.nc+tok.nv+1);
c.ic = c.x(1:tok.nc,:);
c.iv = c.x(tok.nc + (1:tok.nv), :);
c.ip = c.x(tok.nc+tok.nv+1,:);
c.cpasma = c.ip;

% build flux responses
dz = mean(diff(tok.zg));
dr = mean(diff(tok.rg));
[dpsidx_r, dpsidx_z] = gradient(dpsidx, dr, dz);

c.psizr    = dpsidx;
c.psibry   = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rbdef, eqx.zbdef);
c.psicp    = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rcp, eqx.zcp);

c.psixlo   = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rxlo, eqx.zxlo);
c.psixlo_r = gridresponse2pt(tok.rg, tok.zg, dpsidx_r, eqx.rxlo, eqx.zxlo);
c.psixlo_z = gridresponse2pt(tok.rg, tok.zg, dpsidx_z, eqx.rxlo, eqx.zxlo);

c.psixup   = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rxup, eqx.zxup);
c.psixup_r = gridresponse2pt(tok.rg, tok.zg, dpsidx_r, eqx.rxup, eqx.zxup);
c.psixup_z = gridresponse2pt(tok.rg, tok.zg, dpsidx_z, eqx.rxup, eqx.zxup);

c.psixtarglo   = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rxtarglo, eqx.zxtarglo);
c.psixtarglo_r = gridresponse2pt(tok.rg, tok.zg, dpsidx_r, eqx.rxtarglo, eqx.zxtarglo);
c.psixtarglo_z = gridresponse2pt(tok.rg, tok.zg, dpsidx_z, eqx.rxtarglo, eqx.zxtarglo);

c.psixtargup   = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rxtargup, eqx.zxtargup);
c.psixtargup_r = gridresponse2pt(tok.rg, tok.zg, dpsidx_r, eqx.rxtargup, eqx.zxtargup);
c.psixtargup_z = gridresponse2pt(tok.rg, tok.zg, dpsidx_z, eqx.rxtargup, eqx.zxtargup);

c.psitouch     = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rtouch, eqx.ztouch);
c.psitargtouch = gridresponse2pt(tok.rg, tok.zg, dpsidx, eqx.rtargtouch, eqx.ztargtouch);

c.diff_psicp_psibry     = c.psicp - c.psibry;
c.diff_psicp_psitouch   = c.psicp - c.psitouch;
c.diff_psicp_psixlo     = c.psicp - c.psixlo;
c.diff_psicp_psixup     = c.psicp - c.psixup;
c.diff_psicp_psixtarglo = c.psicp - c.psixtarglo;
c.diff_psicp_psixtargup = c.psicp - c.psixtargup;

cdata = c;
































