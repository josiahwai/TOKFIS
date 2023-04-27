function plot_response(tok, x, levels)

if ~exist('levels', 'var'), levels = 10; end

x = x(:);

if length(x) == tok.nc

  r.psi = tok.mpc*x;
  r.Br = tok.gbr2c*x;
  r.Bz = tok.gbz2c*x;
  r.Bp = sqrt(r.Br.^2 + r.Bz.^2);

elseif length(x) == tok.nc+tok.nv

  r.psi = [tok.mpc tok.mpv] * x;
  r.Br = [tok.gbr2c tok.gbr2v] * x;
  r.Bz = [tok.gbz2c tok.gbz2v] * x;
  r.Bp = sqrt(r.Br.^2 + r.Bz.^2);

end

in = inpolygon(tok.rgg, tok.zgg, tok.limdata(2,:), tok.limdata(1,:));
in = in | circshift(in,1,1) | circshift(in,-1,1) | circshift(in,1,2) | circshift(in,-1,2);


figure
for dum = {'psi', 'Br', 'Bz', 'Bp'}

  fd = dum{:};
  data = reshape(r.(fd), tok.nz, tok.nr);
  data(~in) = nan;

  tab = uitab();
  axes('Parent', tab)

  hold on
  contourf(tok.rg, tok.zg, data, levels, 'linewidth', 0.5, 'color', [1 1 1] * 0.8);
  plot_lim(tok, 'k', 'linewidth', 2)
  plot_coils(tok)
  colorbar
  title(fd, 'fontsize', 14)
  
end





% 
%   mpc   = tok.mpc  ;
%   mpv   = tok.mpv  ;
%   mcc   = tok.mcc  ;
%   mvv   = tok.mvv  ;
%   mcv   = tok.mcv  ;
%   gbr2c = tok.gbr2c;
%   gbr2v = tok.gbr2v;
%   gbz2c = tok.gbz2c;
%   gbz2v = tok.gbz2v;
%   nr = tok.nr;
%   nz = tok.nz;
%   rg = tok.rg;
%   zg = tok.zg;
% 
%   mps = [ mpc mpv ];
%   mss = [ mcc mcv ; mcv' mvv];
%   gbr2s = [gbr2c gbr2v];
%   gbz2s = [gbz2c gbz2v];
% 
% 
%   psi = mps*currents;
%   Br_tmp = gbr2s*currents;
%   Bz_tmp = gbz2s*currents;
%   flux = reshape(psi,nr,nz);
%   Br = reshape(Br_tmp,nr,nz);
%   Bz = reshape(Bz_tmp,nr,nz);
%   B = sqrt(Br.^2+Bz.^2);
% 
%   plot_options = struct('iblackbg', 0, 'idxvv', idxvv);
%   plot_sparc_geo(tok, plot_options)
%   hold on
%   [c,h]=contour(rg,zg,flux,clevels);
% %   axis([1,2.9,-1.85,1.85])
%   %if pause_it wait; end
% 
%   % stuff to make the figure tighter
%   % xlabel('')
%   % ylabel('')
%   % ax = gca;
%   % outerpos = ax. OuterPosition;
%   % ti = ax.TightInset;
%   % left = outerpos(1) + ti(1);
%   % bottom = outerpos(2) + ti(2);
%   % ax_width = outerpos(3) - ti(1) - ti(3);
%   % ax_height = outerpos(4) - ti(2) - ti(4);
%   % ax.Position = [left bottom ax_width ax_height];
