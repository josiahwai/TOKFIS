% Some settings for easily loading <spec, init, config> for use with 
% gsdesign. This just populates it with some reasonable values, the user
% will likely want to modify the outputs. 
% 
% INPUTS:
%    shapeclass - string, must be 'LSN', 'USN', or 'LIM' specifying whether
%       a lower single null, upper single null, or limited plasma
%
% OUTPUTS: 
%    spec, init, config for gsdesign
%    

function [spec, init, config] = gsdesign_sparc_easyload(shapeclass, ...
  tok, circ)


% load and define init and shapes
if strcmp(shapeclass, 'LIM')

  init = load('eqLIM').eq;  
  r = init.rbbbs;
  z = init.zbbbs;

  i = r==0 & z==0;      % remove any zero padding
  r(i) = [];
  z(i) = [];

  rx = nan;
  zx = nan;

  [rbdef, i] = min(r);  % specify the touch point as the bdef point
  zbdef = z(i);

elseif strcmp(shapeclass, 'LSN')

  init = load('eqLSN').eq;
  r = init.rbbbs;
  z = init.zbbbs;

  i = r==0 & z==0;      % remove any zero padding
  r(i) = [];
  z(i) = [];

  [zx,i] = min(z);
  rx = r(i);

  rbdef = rx;          % specify the x-point as the bdef point
  zbdef = zx; 

elseif strcmp(shapeclass, 'USN')

  init = load('eqLSN').eq;
  init.psizr = flipud(init.psizr);    % flip everthing vertically to
  init.pcurrt = flipud(init.pcurrt);  %  make it USN instead of LSN
  r = init.rbbbs;
  z = -init.zbbbs;

  i = r==0 & z==0;      % remove any zero padding
  r(i) = [];
  z(i) = [];
  
  [zx,i] = max(z);
  rx = r(i);

  rbdef = rx;           % specify the x-point as the bdef point
  zbdef = zx;

end

% remove fields from init that can cause issues with gsdesign
fds = {'cc', 'cc2', 'ic', 'icx', 'iv', 'ivx'};
for i = 1:length(fds)
  x = fds{i};
  if isfield(init, x)
    init = rmfield(init, x);
  end
end


% sort and interpolate boundary
[r,z] = sort_ccw(r,z);
r(end+1) = r(1);
z(end+1) = z(1);
[r,z] = interparc(r,z,100,0,0);

% load tok into config
config = tok;   

% constrain profiles to three degrees of freedom
config.constraints = 1;  

% a generous number of knots ensures profiles can be made
config.psikn = (0:config.nr-1)/(config.nr-1);  


% specify boundary
spec.targets.rsep = r;
spec.targets.zsep = z;
spec.weights.sep = ones(size(spec.targets.rsep)) * 1e3;

% specify plasma current
spec.targets.cpasma = 3e6;
spec.weights.cpasma = 1e-2;

% specify boundary flux
spec.targets.psibry = 0;
spec.weights.psibry = 10;

% target zero currents
spec.targets.ic = zeros(tok.nc,1);
spec.weights.ic = ones(size(spec.targets.ic)) * 1e-5;

% specify li
spec.targets.li = 1;
spec.weights.li = 10;

% specify betap
spec.targets.betap = 0.4;
spec.weights.betap = 1e2;


% specify x-points
if ~isnan(rx)
  spec.targets.rx = rx;
  spec.targets.zx = zx;
  spec.weights.x  = 100;
end

% specify boundary-defining point
spec.targets.rbdef = rbdef;
spec.targets.zbdef = zbdef;
spec.weights.bdef = 1;


% specify circuit limits
spec.limits.ic = [circ.icmin(:) circ.icmax(:)];

% a reasonable number of iterations
spec.mxiter = 20;    

















