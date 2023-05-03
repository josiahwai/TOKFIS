
function u = shape_cinv_controller(data, datadot, tnow, CONFIG)

persistent ei tprev
if isempty(ei), ei = data*0; end
if isempty(tprev), tprev = tnow; end

dt = 1e-4;
r = interp1(CONFIG.r.Time, CONFIG.r.Data, tnow, 'linear', 'extrap');
rnext = interp1(CONFIG.r.Time, CONFIG.r.Data, tnow+dt, 'linear', 'extrap');
r = r(:);
rnext = rnext(:);
rdot = (rnext - r)/dt;


e = r-data;                % errors
ei = ei + e*(tnow-tprev);  % integral of errors
edot = rdot - datadot;
u = zeros(8,1);
ip = CONFIG.ip;
iy = CONFIG.iy;


% radial control
vec = [0 0 0 0 1 0 0 0]';
kp = 8e-3;
ki = 1e-2;
kd = 2e-4;
c = 0.5;
tmp = kp*e(iy.rcur) + ki*ei(iy.rcur) + kd * (c*rdot(iy.rcur) - datadot(iy.rcur));
u = u + vec * ip * tmp;

 
% vertical control
vec = [0 0.2 0.2 0.8 0 -0.8 -0.2 -0.2]';
kp = 0.0021;
ki = 0;
kd = 1.2e-4;
tmp = kp*e(iy.zcur) + ki*ei(iy.zcur) + kd*edot(iy.zcur);
u = u + vec * ip * tmp;


% Ip control
vec = [1 0 0 0 0 0 0 0]';
kp = -0.1;
ki = -0.4;
kd = -0.001;
c = 0.5;
tmp = kp*e(iy.cpasma) + ki*ei(iy.cpasma) + kd * (c*rdot(iy.cpasma) - datadot(iy.cpasma));
u = u + vec * tmp;

tprev = tnow;



% % SHAPE CONTROL
% fds = {'psixlo', 'diff_psicp_psixlo', 'psixtarglo_r','psixtarglo_z'};
% 
% targ = CONFIG.targ;
% 
% cdata = build_cdata(tok.mpc, tok, targ);
% C = struct2vec(cdata, fds);
% 

% eq = = measure_y(y, targ, CONFIG);
% 
% 
% 
% x = zeros(45,1);
% cv = CONFIG_DATA.cv;
% % fds = fieldnames(cv.iy);
% tok = CONFIG_DATA.tok;
% targ = CONFIG_DATA.targ;
% c = CONFIG_DATA.c;
% circ = CONFIG_DATA.circ;
% 
% dpsizrdx = squeeze(interp1(CONFIG_DATA.dpsizrdx.Time, CONFIG_DATA.dpsizrdx.Data, tnow, 'linear', 'extrap'));
% dpsizrdx = dpsizrdx(:,circ.iic);
% 
% eq.psizr = reshape(data(c.iy.psizr), 33, 33);
% ref.rcp        = interp1(targ.rcp.Time, targ.rcp.Data, tnow, 'linear', 'extrap');
% ref.zcp        = interp1(targ.zcp.Time, targ.zcp.Data, tnow, 'linear', 'extrap');
% ref.rxtargup   = interp1(targ.rxup.Time, targ.rxup.Data, tnow, 'linear', 'extrap');
% ref.zxtargup   = interp1(targ.zxup.Time, targ.zxup.Data, tnow, 'linear', 'extrap');
% ref.rxtarglo   = interp1(targ.rxlo.Time, targ.rxlo.Data, tnow, 'linear', 'extrap');
% ref.zxtarglo   = interp1(targ.zxlo.Time, targ.zxlo.Data, tnow, 'linear', 'extrap');
% ref.rtargtouch = min(tok.limdata(2,:));
% ref.ztargtouch = 0;
% ref.rbdef = data(c.iy.rbdef);
% ref.zbdef = data(c.iy.zbdef);
% ref.rtouch = 0.315;
% ref.ztouch = 0;
% [ref.rxlo, ref.zxlo] = isoflux_xpFinder(eq.psizr, 0.6, -1, tok.rg, tok.zg);
% [ref.rxup, ref.zxup] = isoflux_xpFinder(eq.psizr, 0.6,  1, tok.rg, tok.zg);
% ref.rbnt = 0.315;
% ref.zbnt = 0;
% 
% if tnow >= 0.22
%   ref.rbnt = ref.rxlo;
%   ref.zbnt = ref.zxlo;
% else
%   ref.rbnt = 0.315;
%   ref.zbnt = 0;
% end
% 
% y = measure_y(data, ref, CONFIG_DATA);
% yk = struct2vec(y, fds);
% yr = structts2vec(CONFIG_DATA.targ, fds, tnow);
% e = yr - yk;
% e(isnan(e)) = 0;
% dx0 = e(cv.iy.icx);
% if isempty(dx1)
%   dx1 = dx0;
% end
% 
% ts = 5e-6;
% dt = 0.001;
% if tnow>=0.01 && floor((tnow+ts/2)/dt)>floor((tnow-ts/2)/dt)
%   
%   cdata = build_cmat(dpsizrdx, ref, circ, tok);
%   C = struct2vec(cdata, fds);
%   C = C(:,circ.iicx);    
% 
%   %% LIMITED
%   
%   wy.ip = 0;
%   wy.icx = [0 0 0 0 0 0 0 0]' * 0;
%   wy.diff_psicp_psitouch = ones(length(ref.rcp),1) * 3e10;
%   wy.diff_psicp_psixlo   = ones(length(ref.rcp),1) * 0;
%   wy.diff_psixlo_psixtarglo = 0;
%   wy.psixtarglo_r = 1e8;
%   wy.psixtarglo_z = 1e8;
%   wx = [1e3 1 1 1 1 1 1 1]';
%   
%   Wy = diag(struct2vec(wy,fds));
%   Wx = diag(wx);
%   
%   if 1
%     Dx = eye(8);
%     Dx(4,6) = 1;
%     Dx(6,4) = 1;
%     Dx(6,6) = -1;
%     wx = [1e3 1 1 1 1 1e2 1 1];    
%     Wx = Dx'*diag(wx)*Dx;
%   end
%           
%   H = C'*Wy*C + Wx;
%   f = -C'*Wy*e;
%   
%   % boundary constraints
%   A = zeros(3,8);
%   b = ones(3,1);
%   
%   s = 0.01; % increasing this pushes it to be more LSN than USN
%   A(1,:) = -cdata.diff_psixlo_psixup(circ.iic);
%   b(1) = y.diff_psixlo_psixup - s;
%   
%   s = 0.005; % increasing this pushes it to be more TOUCH than LSN
%   A(2,:) = cdata.diff_psixlo_psitouch(circ.iic);
%   b(2) = -y.diff_psixlo_psitouch - s;
%   
%   s = 0.005; % increasing this pushes it to be more TOUCH than USN
%   A(3,:) = cdata.diff_psixup_psitouch(circ.iic);
%   b(3) = -y.diff_psixup_psitouch - s;
%    
%   ic0 = data(c.iy.ic);
%   slack = 0; 
%   
%   % restrict PF3U/L to help with VS
%   circ.limits.ic([4 6], 2) = interp1([0 0.3 1], [3000 -2000 -2000], tnow); 
%   ub = circ.limits.ic(:,2) - ic0 + slack;
%   lb = circ.limits.ic(:,1) - ic0 - slack;  
%   %   ub = inf(8,1);
%   %   lb = -inf(8,1);
%   opts = optimoptions('quadprog', 'algorithm', 'active-set', 'Display', 'off');
%   x0 = -H\f;
%   dxl = quadprog(H,f,A,b,[],[],lb,ub,x0,opts);
% 
% 
%   %% DIVERTED
%   
%   wy.diff_psicp_psitouch = ones(length(ref.rcp),1) * 0;
%   wy.diff_psicp_psixlo   = ones(length(ref.rcp),1) * 3e10;
%   % wy.diff_psicp_psixlo(ref.zcp>0.5 | ref.zcp<-0.5) = 3e11;
%   k = [12:14 26:28];
%   wy.diff_psicp_psixlo(k) = 3e11;
%   wy.diff_psixlo_psixtarglo = 0;
%   wy.psixtarglo_r = 6e9;
%   wy.psixtarglo_z = 6e9;
%   wy.diff_psixlo_psixtarglo = 0;
%   
%   Wy = diag(struct2vec(wy,fds));
%   
%   H = C'*Wy*C + Wx;
%   f = -C'*Wy*e;
%   
%   % boundary constraints
%   A = zeros(3,8);
%   b = ones(3,1);
%   
%   s = 0.01; % increasing this pushes it to be more LSN than USN
%   A(1,:) = -cdata.diff_psixlo_psixup(circ.iic);
%   b(1) = y.diff_psixlo_psixup - s;
%   
%   s = 0.02; % increasing this pushes it to be more LSN than TOUCH
%   A(2,:) = -cdata.diff_psixlo_psitouch(circ.iic);
%   b(2) = y.diff_psixlo_psitouch - s;
%   
%   
%   x0 = -H\f;
%   dxd = quadprog(H,f,A,b,[],[],lb,ub,x0,opts); 
%   
%   %% BLEND  
%   alpha = blend(tnow, 0.18, 0.25, 2);
%   dx1 = (1-alpha)*dxl + alpha*dxd;        
% end
% 
% dx1(5) = 0;
% alpha = blend(tnow, 0.05, 0.18, 2);
% dx = (1-alpha)*dx0 + alpha*dx1;
% 
% dxr = dx(2:end);
% dpsiin  = alpha * (y.psibry - y.psicp(21)); 
% dpsiout = alpha * (y.psibry - y.psicp(1));
% 
% if dpsiout < 0
%   dpsiout = dpsiout * 2;
% end
% 
% 
% if 1
%   
%   psizr0 = data(c.iy.psizr);
%   psizr0 = reshape(psizr0, c.nz, c.nr);
%   psibry0 = y.psibnt;
%   dpsibrydx = zeros(1,8);
%   for i = 1:8
%     dpsibrydx(i) = bicubicHermite(c.rg, c.zg, reshape(dpsizrdx(:,i), c.nz, c.nr), ref.rbnt, ref.zbnt);
%   end
%   psizr = psizr0 + reshape(dpsizrdx*dx, c.nz, c.nr);
%   psibry = psibry0 + dpsibrydx * dx;
%   
%   rcp = interp1(targ.rcp.Time, targ.rcp.Data, tnow, 'linear', 'extrap');
%   zcp = interp1(targ.zcp.Time, targ.zcp.Data, tnow, 'linear', 'extrap');
%   rxtarglo   = interp1(targ.rxlo.Time, targ.rxlo.Data, tnow, 'linear', 'extrap');
%   zxtarglo   = interp1(targ.zxlo.Time, targ.zxlo.Data, tnow, 'linear', 'extrap');
%   
%   if 1
%     hold on
%     contour(c.rg, c.zg, psizr0, [psibry0 psibry0], 'b')
%     contour(c.rg, c.zg, psizr, [psibry psibry], '--r')
%     scatter(rcp, zcp, 'r', 'filled')
%     plot(rxtarglo, zxtarglo, 'bx', 'linewidth', 3, 'markersize', 12)
%     axis equal
%     axis([0 1.8 -1.8 1.8])
%     set(gcf,'Position', [116 178 337 491]); 
%     psixlo  = bicubicHermite(c.rg, c.zg, psizr, ref.rxlo, ref.zxlo);
%     psixup  = bicubicHermite(c.rg, c.zg, psizr, ref.rxup, ref.zxup);
%     psitouch  = bicubicHermite(c.rg, c.zg, psizr, 0.315, 0);
%   end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 























