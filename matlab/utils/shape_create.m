% PURPOSE: generate an equilibrium shape from low-order shape parameters.
%
% INPUTS:
%   npts - size of output (r,z) vectors
%   s - struct with geometric parameters such as elongation and
%       triangularity. s may contain any or all of the following
%       (any non-specified parameters are overwritten by defaults):
%
%  rsurf  - r of geometric center
%  zsurf  - z of geometric center
%  aminor - plasma radius
%  elong  - plasma elongation
%  triu   - upper triangularity
%  tril   - lower triangularity
%  squo   - upper outer squareness
%  sqlo   - lower outer squareness
%  squi   - upper inner squareness
%  sqli   - lower inner squareness
%  c_xplo - Parameter that affects lower x-point position wrt main plasma
%           and the x-point angle (reasonable range is 0-0.5). If <=0,
%           no x-point is formed.
%  c_xpup - same for upper x-point
%
%
%  EXAMPLE:
%  s.triu = 0.5;
%  npts = 200;
%  [r,z,s] = shape_create(s, npts)
%  plot(r,z)
%  axis equal
%
%  NOTE: specifying x-points and the shaping parameters together is usually
%  not self-consistent, which means the output 's' parameters will only 
%  approximately match the input 's' if x-points are specified.  

function [r,z,s] = shape_create(s, npts)

default.rsurf = 1.1;
default.zsurf = 0;
default.aminor = 0.5;
default.elong = 1.8;
default.triu = 0.3;
default.tril = 0.4;
default.squo = -0.1;
default.sqlo = -0.2;
default.squi = -0.1;
default.sqli = -0.3;
default.c_xplo = 0.1;
default.c_xpup = 0;

s = copyfields(default, s, [], 1);

% make a circle
th = linspace(0, 2*pi, 400);
x = cos(th);
y = sin(th);

% add the x-points
if s.c_xplo >= 0
  x(end+1) = 0;
  y(end+1) = min(y) - s.c_xplo;
end

if s.c_xpup >= 0
  x(end+1) = 0;
  y(end+1) = max(y) + s.c_xpup;
end

% create a convex hull from the circle + x-points
i = convhull(x,y);
x = x(i);
y = y(i);
[x,y] = interparc(x,y,500,0,0);  % generate higher point density

% apply shaping parameters
[r,z] = shape_edit(x,y,s);
s = shape_analysis(r,z);

% interpolate and sort
[r,z] = interparc(r,z,npts,0,0);
[r,z] = sort_ccw(r,z);













































