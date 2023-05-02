function rzip_analysis(models, Kr, Kz, Kip, delay, tok, opt)

% set defaults
default.rstep = 0.01;
default.zstep = 0.01;
default.ipstep = 1e4;
default.t = linspace(0, 0.3, 500);
if ~exist('opt','var'), opt = struct; end
opt = copyfields(default, opt, [], 1);



% close loops with RZIP controllers
s = tf('s');
ssdelay = pade(exp(-delay*s), 4);
P = {};
ip = {};
Pcl = {};

for i = 1:length(models)

  sys = models{i};
  ip = sys.eq.cpasma;

  I = eye(size(sys.amat));
  Cip = I(end,:);
  Cz = sys.dzcurdx;
  Cr = sys.drcurdx;

  Crzip = [Cr; Cz; Cip];

  P{i} = ss(sys.amat, sys.bmat, I, 0);

  Krzip = [Kr*ip Kz*ip Kip];


  L{i} = Crzip * P{i} * ssdelay * Krzip;
  Pcl{i} = feedback(L{i}, eye(3));  
  Gcl{i} = feedback(P{i} * Krzip, Crzip * ssdelay);

end


%% Everything below is plotting
t = opt.t;

co = mat2cell(colororder, ones(7,1), 3);   % helps with defining colors

labels = {'Rcur', 'Zcur', 'Ip'};

y = [];
x = [];

for i = 1:length(models)  
  stepsize = [opt.rstep opt.zstep opt.ipstep];
  y(i,:,:,:) = step(Pcl{i} * diag(stepsize), t);
  x(i,:,:,:) = step(Gcl{i} * diag(stepsize), t);
end
ic = x(:,:,1:tok.nc,:);


% make figure
fig = figure;
fig.Position = [523 109 805 851];
co = mat2cell(colororder, ones(7,1), 3);   % helps with defining colors

for i = 1:3

  tab = uitab();
  axes('Parent', tab);
  sgtitle([labels{i} ' response'], 'fontsize', 20)
  hold on

  % Plot RZIp
  for j = 1:3
    subplot(2,2,j)
    hold on
    if i == j
      plot(t, y(:,:,j,i), 'color', co{1}, 'linewidth', 1.5)
    else
      plot(t, y(:,:,j,i), 'color', co{2}, 'linewidth', 0.5)
    end
    str = sprintf('From %s to %s', labels{i}, labels{j});
    title(str, 'fontsize', 14)
  end

   
  % Nyquist plot of loop gain
  p = uipanel('Parent', tab, 'Position', [0.5 0.05 0.45 0.45], 'units', 'normalized', 'bordertype', 'none');
  axes('Parent', p)
  hold on
  w = {1e-1, 1e4};
  for k = 1:length(models)
    nyquist(L{k}(i,i), w, 'b')
    axis([-4 1 -2 2])
    drawnow
  end
  title(['Nyquist ' labels{i}])

  

  % Plot coil currents
  tab = uitab();
  axes('Parent', tab);
  hold on
  tilelayout = tiledlayout(tab, 'flow');

  for j = 1:tok.nc
    ax = nexttile(tilelayout);
    hold on
    plot(ax, t, ic(:,:,j,i), 'color', co{1})
    title(tok.ccnames{j}, 'fontsize', 12)
    xlabel('Time [s]', 'fontsize', 12)
    ylabel('Current [A]', 'fontsize', 12)
  end
  sgtitle([labels{i} ' coil responses'], 'fontsize', 20)
  
end


































