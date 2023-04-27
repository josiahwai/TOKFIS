
function m = gridresponse2pt(rg, zg, x, r, z)

nx = size(x,2);
r = r(:);
z = z(:);
ny = length(r);

m = zeros(ny,nx);

for i = 1:nx
  xi = reshape(x(:,i), length(zg), length(rg));
  m(:,i) = bicubicHermite(rg, zg, xi, r, z);
end











