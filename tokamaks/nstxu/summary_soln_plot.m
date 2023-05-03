% Creates a ui figure for plotting equilibrium solutions

function fig = summary_soln_plot(y, iy, tok, targ)

if ~exist('targ', 'var')
  targ.rcp = nan;
  targ.zcp = nan;
end

% create figure
fig = figure;
fig.Position = [754   322   425   611];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];


% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = y.Time(1);
s.Max = y.Time(end);
s.Value = y.Time(1);
s.Callback = {@sliderCallback, y, iy, tok, targ};

% text edit button
e = uicontrol(fig, 'style', 'edit');
e.Units = 'normalized';
e.Position = [0.8 0.07 0.15 0.05];
e.Callback = {@editCallback, y, iy, tok, targ};

plot_shape(0, y, iy, tok, targ)


% slider callback
function sliderCallback(src, event, y, iy, tok, targ)
  t = src.Value;
  plot_shape(t, y, iy, tok, targ)
end


% text edit callback
function editCallback(src, event, y, iy, tok, targ)
  t = str2double(src.String);
  plot_shape(t, y, iy, tok, targ)
end


% plot shape targets
function plot_shape(t, y, iy, tok, targ)

  [~,i] = min(abs(t-y.Time));
  
  psizr = y.Data(iy.psizr, i);
  psizr = reshape(psizr, tok.nz, tok.nr);
  [rx,zx,psix] = isoflux_xpFinder(tok.rg, tok.zg, psizr, 0.6, -1);

  cla
  hold on  
  plot_lim(tok, 'k', 'linewidth', 2)
  contour(tok.rg, tok.zg, psizr, 10, 'color', [1 1 1] * 0.8)
  contour(tok.rg, tok.zg, psizr, psix*[1 1], 'r', 'linewidth', 2) 
  scatter(targ.rcp, targ.zcp, 60, 'b', 'filled')
  plot_coils(tok, 'k', 0.2)

  text(-0.25, -0.1, 'Drag slider to view equilibria', 'units', 'normalized', 'fontsize', 11)  
  text(1.05, -0.1, 'Enter time:', 'units', 'normalized', 'fontsize', 11)  
  str = sprintf('Time=%.3f', t);
  title(str, 'fontsize', 14)
  drawnow
  

end


end





































