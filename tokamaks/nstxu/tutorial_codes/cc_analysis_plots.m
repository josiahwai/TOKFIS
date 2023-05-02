function fig = cc_analysis_plots(t, y, u, L, tok, coils2plot)

idx2plot = ccnames2idx(coils2plot, tok);


% make figure
fig = figure;
fig.Position = [523 109 805 851];
co = mat2cell(colororder, ones(7,1), 3);   % helps with defining colors


% plot all coils and voltages
tab = uitab();
axes('Parent', tab);

subplot(211)
title('Coil current step response [A]', 'fontsize', 18)
hold on
grid on
for i = 1:tok.nc
  plot(t, y(:,i,i))
end

subplot(212)
title('Coil voltage step response [V]', 'fontsize', 18)
hold on
grid on
for i = 1:tok.nc
  plot(t, u(:,i,i))
end



% plot responses of individual coils
for i = idx2plot(:)'


  % define tabs and panels to organize figure
  tab = uitab();  
  axes('Parent', tab, 'Visible', 'off') 
  s = sgtitle([tok.ccnames{i} ' Current Response'], 'fontweight', 'bold');

  h = 0.55;
  p1 = uipanel('Parent', tab, 'Position', [0 h 0.5 (1-h-0.07)], 'units', 'normalized', 'bordertype', 'none');
  p2 = uipanel('Parent', tab, 'Position', [0.5 h 0.5 (1-h-0.07)], 'units', 'normalized', 'bordertype', 'none');
  p3 = uipanel('Parent', tab, 'Position', [0 0 1 h], 'units', 'normalized', 'bordertype', 'none');


  % step response of currents
  axes('Parent', p1)
  hold on  
  plot(t, y(:,1:tok.nc,i), 'color', co{2}, 'linewidth', 1)
  plot(t, y(:,i,i), 'color', co{1}, 'linewidth', 2)
  title('Current [A]', 'fontsize', 14)
  ylabel('A')
  xlabel('t')
  
  
  % bode plot of loop gain
  axes('Parent', p2)  
  [mag,phase,w] = bode(L(i,i), {1e-1, 1e4});
  mag = squeeze(mag);
  db = mag2db(mag);
  phase = squeeze(phase);  
  subplot(211)
  semilogx(w, db)
  grid on
  hold on
  str = sprintf('Bode diagram: %s loop gain', tok.ccnames{i});
  title(str)
  yline(0)
  xlabel('')
  ylabel('Magnitude (dB)')  
  subplot(212)
  semilogx(w, phase)
  grid on
  hold on
  xlabel('Frequency (rad/s)')
  ylabel('Phase (deg)')


  % plot cross currents
  axes('Parent', p3)
  nrows = 2;
  ncols = ceil(tok.nc/nrows);

  for j = 1:tok.nc
    subplot(nrows, ncols, j)
    hold on
    str = sprintf('From %s to %s', tok.ccnames{i}, tok.ccnames{j});
    title(str, 'fontsize', 12, 'fontweight', 'normal')    
    if i == j
      plot(t, y(:,j,i), 'color', co{1}, 'linewidth', 2)
    else
      plot(t, y(:,j,i), 'color', co{2})     
    end
    ylabel('A')  
  end
 

  % new tab
  tab = uitab();  
  axes('Parent', tab, 'Visible', 'off') 
  s = sgtitle([tok.ccnames{i} ' Voltage Response'], 'fontweight', 'bold');
  h = 0.55;
  p1 = uipanel('Parent', tab, 'Position', [0.25 h 0.5 (1-h-0.07)], 'units', 'normalized', 'bordertype', 'none');  
  p2 = uipanel('Parent', tab, 'Position', [0 0 1 h], 'units', 'normalized', 'bordertype', 'none');


  % step response of voltages
  axes('Parent', p1)
  hold on  
  plot(t, u(:,:,i), 'color', co{2}, 'linewidth', 1)
  plot(t, u(:,i,i), 'color', co{1}, 'linewidth', 2)
  title('Voltage [V]', 'fontsize', 14)
  ylabel('V')
  xlabel('t')


  % plot cross voltages
  axes('Parent', p2)
  for j = 1:tok.nc
    subplot(nrows, ncols, j)
    hold on
    str = sprintf('From %s to %s', tok.ccnames{i}, tok.ccnames{j});
    title(str, 'fontsize', 12, 'fontweight', 'normal')    
    if i == j
      plot(t, u(:,j,i), 'color', co{1}, 'linewidth', 2)
    else
      plot(t, u(:,j,i), 'color', co{2})     
    end
    ylabel('V')  
  end

end
