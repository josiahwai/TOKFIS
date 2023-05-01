% obtains plasma response models from several different plasma-response
% codes (gsupdate.m, rzrigid.m, gspert.m, rig.m) and puts them all into a
% common format. 
%
% INPUTS: 
%     eq           - equilibrium 
%     tok          - tok_data_struct geometry
%     iplcirc      - if 1, includes plasma current as a state vector
%     Rp           - plasma resistance in Ohms
%     plasma_model - string that is one of {'gspert', 'rzrig', 'gsupdate', 
%                    'rig', 'vacuum', or 'all'} or a cellstring with any 
%                    combination of these
%
% Restrictions: tested against equilibria created by gsdesign.m, to fit
%    another equilibrium format (e.g. geqdsk file) with gsdesign, recommend
%    using [~,~,eq] = gsdesign_recreate_init(init, tok)
%
% The {'amat', 'bmat', 'dpcurrtdx', 'drcurdx', 'dzcurdx'}
%  

function models = response_models(eq, tok, iplcirc, Rp, plasma_model)

  if ischar(plasma_model)
    plasma_model = cellstr(plasma_model); 
  end
  if any(strcmp(plasma_model, 'all')) 
    plasma_model = {'rzrig', 'gspert', 'gsupdate', 'rig', 'vacuum'};
  end

  N = length(plasma_model);
  models = cell(N,1);

  for i = 1:N
    x = plasma_model{i};  
    % try
      fprintf('building %s sys\n', x);
      fun = str2func(['build_' x '_sys']); 
      models{i} = fun(eq, tok, iplcirc, Rp);
      models{i}.eq = eq;
      models{i}.plasma_model = x;      
    % catch
    % warning(['Could not run ' x ])
    % end  
  end

  if N == 1
    models = models{1};
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is infile functions for calling builds of 
%    gspert/rzrig/gsupdate/rig/vacuum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
       
% =================
% build_gspert_sys
% =================
function sys = build_gspert_sys(eq, tok, iplcirc, Rp)

  b = struct;
  b.tokamak = '';
  b.vacuum_objs = tok;
  b.cccirc = 1:tok.nc;
  b.ichooseq = 4;
  b.irzresp_dynamic = 5;
  b.irzresp_output = 5;
  b.iplcirc = iplcirc;
  b.Rp = Rp;
  b.equil_data = eq;
  s = build_tokamak_system(b);  % builds the dynamics 
  
  % extract and relabel a subset of response derivatives
  r = s.gspert_data;
  if iplcirc
    s.dpcurrtdx = [r.dcphidis r.dcphidip(:)];
    s.drcurdx = [r.drdis r.drdip];
    s.dzcurdx = [r.dzdis r.dzdip];
  else
    s.dpcurrtdx = r.dcphidis;
    s.drcurdx = r.drdis;
    s.dzcurdx = r.dzdis;
  end
  
  fds2copy = {'amat', 'bmat', 'dpcurrtdx', 'drcurdx', 'dzcurdx'};
  sys = copyfields(struct, s, fds2copy, 1);
end




% ===============
% build_rzrig_sys
% ===============
function sys = build_rzrig_sys(eq, tok, iplcirc, Rp)

  b = struct;
  b.tokamak = '';
  b.vacuum_objs = tok;
  b.cccirc = 1:tok.nc;
  b.ichooseq = 4;
  b.irzresp_dynamic = 3;
  b.irzresp_output = 3;
  b.iplcirc = iplcirc;
  b.Rp = Rp;
  b.equil_data = eq;
  s = build_tokamak_system(b);  % builds the dynamics 
  
  % extract and relabel a subset of response derivatives
  r = s.rzrig_data;
  if iplcirc
    s.drcurdx = [r.drdis r.drdip];
    s.dzcurdx = [r.dzdis r.dzdip];  
  else
    s.drcurdx = r.drdis;
    s.dzcurdx = r.dzdis;
  end
  s.dpcurrtdx = r.dcdr(:)*s.drcurdx + r.dcdz(:)*s.dzcurdx;
    
  fds2copy = {'amat', 'bmat', 'dpcurrtdx', 'drcurdx', 'dzcurdx'};
  sys = copyfields(struct, s, fds2copy, 1);
end



% ===================
% build_gsupdate_sys
% ===================
function sys = build_gsupdate_sys(eq, tok, iplcirc, Rp)

  c = gsconfig(tok);
  c.evolve_option = 1;
  [d,~,r] = gsupdate(eq.x, eq, c, 1);
  sys = struct;

  dvdx = [eye(c.nic+c.niv, c.nic+c.niv+3); r.dcpasmadx; r.dlidx; r.dbetapdx];
  dxdv = inv(dvdx); 
  dpcurrtdx = r.dpcurrtdx * dxdv;
  dpcurrtdx(:,end-2:end) = [];

  X = [tok.mpc tok.mpv]' * dpcurrtdx;
  M = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
  R = diag([tok.resc; tok.resv]);
  lstar = M+X;
  lstari = inv(M+X);
  amat = -lstari * R;
  bmat = lstari(:,1:tok.nc);

  
  % plasma current evolution from gsupdate is weird due to the usage of 
  % sf and sp coefficients as the states, lets put in something more reasonable:
  
  if iplcirc
    lstari = -amat * inv(R);
    lstar = inv(lstari);
    [Lp, ~, ~, ~, mcIp, mvIp] = inductance(eq, tok);
    lstar = [lstar [mcIp; mvIp]; [mcIp' mvIp' Lp]];
    R = blkdiag(R, Rp);
    lstari = inv(lstar);
    amat = -lstari * R;
    bmat = lstari(:,1:tok.nc);
  end
  
  sys.amat = amat;
  sys.bmat = bmat;

  % map responses from x:=[ic; iv; cpasma; sf; sp] to v:=[ic; iv; cpasma; li; betap] 
  dvdx = [eye(c.nic+c.niv, c.nic+c.niv+3); r.dcpasmadx; r.dlidx; r.dbetapdx];
  dxdv = inv(dvdx); 
  if iplcirc
    dxdv = dxdv(:,1:end-2);
  else
    dxdv = dxdv(:,1:end-3);
  end
  sys.dpcurrtdx = r.dpcurrtdx * dxdv;
  sys.drcurdx   = r.drcurdx * dxdv;
  sys.dzcurdx   = r.dzcurdx * dxdv;
  
end



% =============
% build_rig_sys
% =============
function sys = build_rig_sys(eq, tok, iplcirc, Rp)

  r = rig(eq, tok);
  
  sys = struct;
  
  if iplcirc
    [Lp, ~, ~, ~, mcIp, mvIp] = inductance(eq, tok);
    lstar = [r.lstar [mcIp; mvIp]; [mcIp' mvIp' Lp]];
    R = diag([r.res; Rp]);
    lstari = inv(lstar);  
    sys.amat = -lstari*R;
    sys.bmat = lstari(:,1:tok.nc);
  
    dpcurrtdip = eq.pcurrt(:) / eq.cpasma;   % assume Ip response is a uniform scaling
    sys.dpcurrtdx = [r.djdx dpcurrtdip];
    sys.drcurdx = [r.drcdx 0];
    sys.dzcurdx = [r.dzcdx 0];
  
  else
    sys.amat = r.amat;
    sys.bmat = r.bmat;
    sys.dpcurrtdx = r.djdx;
    sys.drcurdx = r.drcdx;
    sys.dzcurdx = r.dzcdx;
  end
end


% ================
% build_vacuum_sys
% ================
function sys = build_vacuum_sys(eq, tok, iplcirc, Rp)


  lstar = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
  R = diag([tok.resc; tok.resv]);
  
  if iplcirc
    [Lp, ~, ~, ~, mcIp, mvIp] = inductance(eq, tok);
    lstar = [lstar [mcIp; mvIp]; [mcIp' mvIp' Lp]];
    R = blkdiag(R, Rp);
  end
  lstari = inv(lstar);
  nx = size(lstar,1);
  
  sys = struct;
  sys.amat = -lstari*R;
  sys.bmat = lstari(:,1:tok.nc);
  sys.dpcurrtdx = zeros(tok.nz*tok.nr, nx);
  sys.drcurdx = zeros(1,nx);
  sys.dzcurdx = zeros(1,nx);

end

















































