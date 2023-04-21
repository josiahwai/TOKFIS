% sorts the points in (r,z) according to increasing counterclockwise angle
% that they make with (r0,z0). 

function [r,z] = sort_ccw(r,z,r0,z0)

if nargin == 2
  warning('off', 'MATLAB:polyshape:repairedBySimplify')
  P = polyshape(r,z);
  warning('on', 'MATLAB:polyshape:repairedBySimplify')
  [r0,z0] = centroid(P);  
end

th = atan2(z-z0, r-r0);
th = wrapTo2Pi(th);
[~,i] = sort(th);
r = r(i);
z = z(i);







