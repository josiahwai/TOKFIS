function vs_analysis(models, Kz, delay, opt)

stepsize = 0.01;
t = linspace(0, 0.3, 500);

s = tf('s');
ssdelay = pade(exp(-delay*s), 4);

N = length(models);
CL = cell(N,1);
L = cell(N,1);

ic = [];
zcur = [];

for i = 1:N

  sys = models{i};
  ip = sys.eq.cpasma;
  G = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);
  
  if opt.reduce_model
    G = reduce(G, opt.model_order);  
  end

  L{i} = sys.dzcurdx * G * ssdelay * Kz * ip;       % loop gain
  CL{i} = feedback(L{i}, 1);                        % closed loop from zcur r to y
  CLX{i} = feedback(G*Kz*ip, sys.dzcurdx*ssdelay);  % closed loop gain from zcur to x


  x = step(CLX{i}, t) * stepsize;

  ic(i,:,:) = x(:, 1:tok.nc);
  zcur(i,:) = step(CL{i}, t) * stepsize;
end


%% Everything else is plotting


% make figure
fig = figure;
fig.Position = [523 109 805 851];
co = mat2cell(colororder, ones(7,1), 3);   % helps with defining colors

tab = uitab();

% Step responze of zcur
p = uipanel('Parent', tab, 'Position', [0.2 0.55 0.6 0.45], 'units', 'normalized'); % , 'bordertype', 'none');
axes('Parent', p)
hold on
grid on
title('Step Response', 'fontsize', 16)
plot(t, zcur)
ymax = 2 * quantile(zcur(:), 0.9);
ymin = 0;
slack = (ymax - ymin) * 0.3;
ylim([ymin-slack ymax+slack])
xlabel('Time [s]', 'fontsize', 16)
ylabel('Z cur [m]', 'fontsize', 16)


% Nyquist plot of loop gain
p = uipanel('Parent', tab, 'Position', [0.2 0.05 0.6 0.45], 'units', 'normalized'); % , 'bordertype', 'none');
axes('Parent', p)
hold on
for i = 1:N
  nyquist(L{i})
  axis([-4 1 -2 2])
  drawnow
end


% Plot coil currents
tab = uitab();
tilelayout = tiledlayout(tab, 'flow');
for i = 1:tok.nc
  ax = nexttile(tilelayout);
  hold on
  plot(ax, t, ic(:,:,i))
  title(tok.ccnames{i}, 'fontsize', 16)
  xlabel('Time [s]', 'fontsize', 16)
  ylabel('Current [A]', 'fontsize', 16)
end
  






























