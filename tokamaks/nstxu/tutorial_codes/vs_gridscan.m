
function [gamma, kps, kds] = vs_gridscan(sys, vec, kp0, kd0, tau, opt)



% set default options
default.nkp = 15;
default.nkd = 15;
default.mag_up = 2;
default.mag_lo = -2;
default.plotit = 1;
default.model_order = 40;
default.pade_order = 4;
if ~exist('opt','var'), opt = struct; end
opt = copyfields(default, opt, [], 1);


% define grid to search over
kps = kp0 * logspace(opt.mag_lo, opt.mag_up, opt.nkp);
kds = kd0 * logspace(opt.mag_lo, opt.mag_up, opt.nkd);
gamma = nan(opt.nkd, opt.nkp);


% dynamics stuff
A = sys.amat;
Bz = sys.bmat * vec;
Cz = sys.dzcurdx;
I = eye(size(A));
ip = sys.eq.cpasma;
s = tf('s');
P = ss(A, Bz, Cz*ip, 0);   % open loop plant
P = reduce(P, opt.model_order); 




% search grid for stable values
for i = 1:length(kps)
  for j = 1:length(kds)
    kp = kps(i);
    kd = kds(j);

    K = (kp + kd*s) * exp(-tau*s);  % PD controller with time delay
    K = pade(K, opt.pade_order);
    CL = feedback(P*K, 1);
    gamma(j,i) = max(real(pole(CL))); % max(real(pole(pade(CL,2))));

  end
end


x = log10(kps);
y = log10(kds);
[x,y] = meshgrid(x,y);

i = find(gamma<0);
j = boundary(x(i), y(i));
k = i(j);
xs = x(k);   % xs and ys are the boundaries of the stable (kp,kd)
ys = y(k); 

if isempty(xs)
  fprintf('\n\nwarning vs_gridscan.m: No stable region found.\n\n')
else
  grid off
  hold on
  lim = [min(x(:)) max(x(:)) min(y(:)) max(y(:))];
  co = mat2cell(colororder, ones(7,1), 3);
  h = fillout(xs, ys, lim);
  h.FaceColor = co{1};
  h.FaceAlpha = 0.05;
  plot(xs, ys, 'k', 'linewidth', 0.5)
  axis(lim)
  xlabel('log10 (kp)', 'fontsize', 16)
  ylabel('log10 (kd)', 'fontsize', 16)
  title('Vertical Controller Stable Regions', 'fontsize', 14)
  drawnow
end
















































