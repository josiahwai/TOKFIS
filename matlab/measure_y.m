
function y = measure_y(eq, ref, CONFIG_DATA)





% c = CONFIG_DATA.c;
% psizr = reshape(data(c.iy.psizr), c.nr, c.nz);
% targ = CONFIG_DATA.targ;
% circ = CONFIG_DATA.circ;
% 
% fds1 = {'psixtarglo', 'psixtarglo_r', 'psixtarglo_z', 'psixtargup', ...
%   'psixtargup_r', 'psixtargup_z', 'psixlo', 'psixlo_r', 'psixlo_z', ...
%   'psixup', 'psixup_r', 'psixup_z', 'psitouch', 'psibry', 'diff_psixlo_psixup', ...
%   'diff_psixtarglo_psixtargup', 'psitargtouch', 'psibnt', ...
%   'diff_psibry_psibnt', 'diff_psixlo_psitouch', 'diff_psixup_psitouch', ...
%   'diff_psixlo_psixtarglo', 'diff_psixup_psixtargup'};
% 
% fds2 = {'psicp', 'diff_psicp_psixtarglo', 'diff_psicp_psixtargup', ...
%   'diff_psicp_psitargtouch', 'diff_psicp_psixlo', 'diff_psicp_psixup', ...
%   'diff_psicp_psitouch', 'diff_psicp_psibry', 'diff_psicp_psibnt'};
% 
% dims1 = {1,1};
% dims2 = {length(ref.rcp),1};
% 
% y = struct;
% y = initializeStructDim(y, fds1, dims1);
% y = initializeStructDim(y, fds2, dims2);
% 
% % currents
% y.ip = data(c.iy.cpasma);
% y.icx = data(c.iy.ic);
% y.ivx = CONFIG_DATA.models.bal.Tv * data(c.iy.iv);
% 
% % flux at control points
% rcp = min(max(ref.rcp,c.rg(3)), c.rg(end-2));
% zcp = min(max(ref.zcp,c.zg(3)), c.zg(end-2));
% y.psicp = bicubicHermite(c.rg, c.zg, psizr, rcp(:), zcp(:));
% y.psibry = data(c.iy.psibry);
% 
% % flux at target touch/x-points
% [y.psixtarglo, y.psixtarglo_r, y.psixtarglo_z] = bicubicHermite(c.rg, c.zg, psizr, ref.rxtarglo, ref.zxtarglo);
% [y.psixtargup, y.psixtargup_r, y.psixtargup_z] = bicubicHermite(c.rg, c.zg, psizr, ref.rxtargup, ref.zxtargup);
% [y.psixlo, y.psixlo_r, y.psixlo_z] = bicubicHermite(c.rg, c.zg, psizr, ref.rxlo, ref.zxlo);
% [y.psixup, y.psixup_r, y.psixup_z] = bicubicHermite(c.rg, c.zg, psizr, ref.rxup, ref.zxup);
% y.psitouch = bicubicHermite(c.rg, c.zg, psizr, ref.rtouch, ref.ztouch);
% % y.psibnt = bicubicHermite(c.rg, c.zg, psizr, ref.rbnt, ref.zbnt);
% 
% 
% 
% % control points vs target touch/x-points
% ONE = ones(length(rcp), 1);
% 
% y.diff_psicp_psixtarglo = y.psicp - ONE * y.psixtarglo;
% y.diff_psicp_psixtargup = y.psicp - ONE * y.psixtargup;
% y.diff_psicp_psibry = y.psicp - ONE * y.psibry;
% % y.diff_psicp_psibnt = y.psicp - ONE * y.psibnt;
% y.diff_psicp_psitouch = y.psicp - ONE*y.psitouch;
% y.diff_psicp_psixlo = y.psicp - ONE*y.psixlo;
% 
% 
% y.diff_psixlo_psixtarglo = y.psixlo - y.psixtarglo;
% y.diff_psixup_psixtargup = y.psixup - y.psixtargup;
% % y.diff_psibry_psibnt = y.psibry - y.psibnt;
% y.diff_psixlo_psixup = y.psixlo - y.psixup; 
% y.diff_psixlo_psitouch = y.psixlo - y.psitouch;
% y.diff_psixup_psitouch = y.psixup - y.psitouch;


























