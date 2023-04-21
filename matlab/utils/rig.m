% inputs: eq - equilibrium structure, must have fields ic (coil currents),
%              iv (vessel currents), pcurrt (plasma current distribution),
%              and rbbbs, zbbbs the plasma boundary
%
%         tok - tokamak vacuum geometry containing mutual inductances,
%               greens funs, grid, etc. The grid used in eq and tok 
%               must be the same. The size of the mutual inductance 
%               matrices / greens funs must correspond to the size of eq.ic
%               and eq.iv. Must contain the mutual inductances and magnetic
%               field greens funs (mpc, mpv, gbr2c, gbr2v, ...) 
%
% outputs: r - response of the plasma, see r.info
%
% The rigid plasma response model is based on the description and equations
% given in:  M.L. Walker & D.A. Humphreys, "Valid coordinate systems for 
% linearized plasma shape response models in tokamaks", FST, 2006. 

function r = rig(eq, tok)

% read info from tok
nr = tok.nr;
rg = tok.rg;
rgg = tok.rgg;
nz = tok.nz;
zg = tok.zg;
zgg = tok.zgg;
dr = mean(diff(rg));
dz = mean(diff(zg));
mu0 = pi*4e-7;
gbr = [tok.gbr2c tok.gbr2v];
gbz = [tok.gbz2c tok.gbz2v];
mpx = [tok.mpc tok.mpv];

% vacuum fields
x = [eq.ic; eq.iv];
brvac = reshape(gbr * x, nz, nr);
bzvac = reshape(gbz * x, nz, nr);

% current distribution and gradient
J = eq.pcurrt;
ip = sum(J(:));
rc = sum(J(:).*rgg(:)) / ip;
zc = sum(J(:).*zgg(:)) / ip;
[djdr, djdz] = gradient(J, dr, dz);
djdrc = -djdr;                                             % eqn 5
djdzc = -djdz;                                             % eqn 5

% derivatives w.r.t conductors
dfrdx = 2*pi*(rgg(:) .* J(:))' * gbz;                      % eqn 8
dfzdx = -2*pi*(rgg(:) .* J(:))' * gbr;                     % eqn 11

% derivative of vacuum field forces w.r.t centroid
dfrdrc = -2 * pi * sum( rgg(:) .* djdr(:) .* bzvac(:));      % eqn 9

% The for loops here help improve robustness since dfzdzc can be very
% close to zero and (1/dfzdzc) is used later. This does several small 
% shifts of the Br-field spatially, calculates dfzdzc, and then takes the
% median.
dfzdzc = [];
steps = -2:2;
for i = 1:length(steps)
  br = circshift(brvac, steps(i), 1);
  for j = 1:length(steps)
    br = circshift(br, steps(j), 2);
    dfzdzc(i,j) = 2 * pi * sum( rgg(:) .* djdz(:) .* br(:));  % eqn 10
  end
end
dfzdzc = median(dfzdzc(:));


% derivative of hoop force w.r.t centroid
touchinfo = check_if_limited(eq, tok); 

if touchinfo.HFS_limited
  rhfs = min(tok.limdata(2,:));
  a = rc - rhfs;
  dfrhoopdrc = (mu0*ip^2/2) * (1/rc - 1/a);

elseif touchinfo.LFS_limited  
  rlfs = max(tok.limdata(2,:));
  a = rlfs - rc;
  dfrhoopdrc = (mu0*ip^2/2) * (1/rc + 1/a);

else  % diverted
  dfrhoopdrc = mu0*ip^2/(2*rc);                              % eqn 9
end
dfrdrc = dfrdrc + dfrhoopdrc;

% centroid motion
drcdx = -dfrdx / dfrdrc;                                   % eqn 6
dzcdx = -dfzdx / dfzdzc;                                   % eqn 6 

% plasma current response
djdx = djdrc(:) * drcdx + djdzc(:) * dzcdx;           % similar to eqns 3/4

% additional outputs
xmat = mpx' * djdx;                   % eqns 3 & 4 combined

mmat = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
lstar = xmat + mmat;
lstari = inv(lstar);
res = [tok.resc(:); tok.resv(:)];
amat = -lstari * diag(res);
bmat = lstari(:,1:tok.nc);


% documentation and save outputs
info.djdx        = 'response of plasma current distribution to conducting elements (x)';
info.xmat        = 'the x matrix used in dynamics model (lstar = xmat + mmat)';
info.dzcdx       =  'response of current centroid z to conductors';
info.drcdx       = 'response of current centroid z to conductors';
info.mmat        = 'mutual inductance matrix';
info.lstar       = 'plasma-modified mutual inductance';
info.lstari      = 'inv(lstar)';
info.amat        = 'state-space dynamics A matrix';
info.bmat        = 'state-space dynamics B matrix';
info.res         = 'resistivity vector';
info.touchinfo   = 'info on whether and where plasma touches wall';

r = variables2struct(djdx,xmat,dzcdx,drcdx,amat, ...
  bmat, res, lstar, mmat, lstari, touchinfo);

% sort outputs for readability
fds    = fields(r);
[~,i]  = sort(lower(fds));
fds    = fds(i);
r      = reorderstructure(r, fds{:});
info   = reorderstructure(info, fds{:});
r.info = info; 





















