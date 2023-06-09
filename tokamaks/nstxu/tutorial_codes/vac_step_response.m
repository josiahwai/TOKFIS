function [out, figh] = vac_step_response(...
  sys, K, C, tok, circ, stepsize, t, coils, plotit)


filename = 'coilsim_slx.slx';
G = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);
K = ss(K);
Cc = eye(tok.nc, tok.nc+tok.nv);
vmax = circ.vmax;
vmin = circ.vmin;
tmax = max(t);

n = length(coils);
out = cell(n,1);


for i = 1:n

  icoil = circ.iy.(coils{i});  

  ref = zeros(tok.nc,1);
  ref(icoil) = stepsize;
  
  out{i} = run_slx_model(filename, G, K, Cc, vmax, vmin, tmax, ref);  % run simulink model  
  out{i}.x = resample(out{i}.x, t);
  out{i}.y = out{i}.x * C';
end


%% Plotting
figh = [];
if plotit

  figh = figure;
  for i = 1:n
    tab = uitab();
    axes('Parent', tab)
    plot(out{i}.y)
    hold on
    grid on
    title(coils{i}, 'fontsize', 16)
  end
end




















