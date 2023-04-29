function [iuse, gamma] = sys_eig_check(sys, gamma_max, gamma_min)

if ~exist('gamma_max','var'), gamma_max = inf; end
if ~exist('gamma_min','var'), gamma_min = -inf; end

n = length(sys);
iuse = true(n,1);
gamma = nan(n,1);

for i = 1:n

  e = esort(eig(sys{i}.amat));
  gamma(i) = e(1);

  if real(e(2)) > 0
    iuse(i) = 0;
  end

  if real(e(1)) > gamma_max
    iuse(i) = 0;
  end

  if real(e(1)) < gamma_min
    iuse(i) = 0;
  end   
end



































