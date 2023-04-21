% gets some info on if a plasma boundary is limited, whether high- or 
% low-field side, etc.
%
% Inputs: eq  - equilibrium struct with fields (rbbbs, zbbbs) the plasma
%               boundary
%         tok - tokamak geometry struct with field limdata, the limiter
%               geometry

function touchinfo = check_if_limited(eq, tok)

% check if limited and whether on high or low field side
if size(tok.limdata,1) ~= 2
  tok.limdata = tok.limdata';
end
curvexy = [eq.rbbbs(:) eq.zbbbs(:)];
mapxy = [tok.limdata(2,:)' tok.limdata(1,:)'];

% find the nearest point to wall
[xy, dist] = distance2curve(curvexy, mapxy);
[mindist, i] = min(dist);
rnearest = xy(i,1);   
znearest = xy(i,2); 


% decide whether limited
diverted = 1;
limited = 0;
LFS_limited = 0;
HFS_limited = 0;
rtouch = nan;
ztouch = nan;

if mindist < 0.002    % plasma is within 2mm of wall => limited
  rtouch = rnearest;
  ztouch = znearest;
  diverted = 0;
  limited = 1;

  if rnearest > quantile(eq.rbbbs, 0.8)        % touch point is at large radius
    LFS_limited = 1;
  elseif rnearest < quantile(eq.rbbbs, 0.2)    % touch point is at small radius
    HFS_limited = 1;    
  end
end

touchinfo = variables2struct(diverted, limited, LFS_limited, HFS_limited, ...
  rtouch, ztouch, rnearest, znearest);







