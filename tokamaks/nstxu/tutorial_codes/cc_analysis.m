function cc_analysis(K, sys, tok, stepsize, tmax, coils2plot, circ)  

A = sys.amat;
B = sys.bmat;
Cc = eye(tok.nc, tok.nc+tok.nv);
P = ss(A,B,Cc,0);              % plant model (gain from u to y)
G = ss(A,B,eye(size(A)), 0);   % plant model including vessels

K = ss(K);                % controller model
L = P*K;                  % loop transfer function

vmax = circ.vmax;         % power supply max voltage
vmin = circ.vmin;         % power supply min voltage

n = 400;
t = linspace(0, tmax, n);   % timebase for simulation

u = zeros(n, tok.nc, tok.nc);  % initialize outputs to zero
e = zeros(n, tok.nc, tok.nc);
x = zeros(n, tok.nc+tok.nv, tok.nc);

idx2plot = ccnames2idx(coils2plot, tok);  % which coils to simulate

for i = idx2plot(:)'

  ref = zeros(tok.nc,1);   
  ref(i) = stepsize;        % step input on coil i

  
  filename = 'coilsim_slx.slx';
  out = run_slx_model(filename, G, K, Cc, vmax, vmin, tmax, ref);  % run simulink model
  
  out.e = resample(out.e, t);  % resample outputs
  out.x = resample(out.x, t);
  out.u = resample(out.u, t);

  e(:,:,i) = out.e.Data;       % save array of outputs for all coils
  x(:,:,i) = out.x.Data;
  u(:,:,i) = out.u.Data;
end


cc_analysis_plots(t, x, u, L, tok, coils2plot);  % makes lots of plots






















