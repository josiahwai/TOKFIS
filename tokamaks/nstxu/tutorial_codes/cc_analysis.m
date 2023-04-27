function cc_analysis_vlim(K, sys, tok, stepsize, tmax, coils2plot, circ)  

A = sys.amat;
B = sys.bmat;
Cc = eye(tok.nc, tok.nc+tok.nv);
P = ss(A,B,Cc,0);              % plant model (gain from u to y)
G = ss(A,B,eye(size(A)), 0);   % plant model including vessels

K = ss(K);                % controller model
I = eye(size(P));         % identity matrix
L = P*K;                  % loop transfer function


n = 400;
t = linspace(0, tmax, n);
u = zeros(n, tok.nc, tok.nc);
e = zeros(n, tok.nc, tok.nc);
x = zeros(n, tok.nc+tok.nv, tok.nc);

idx2plot = ccnames2idx(coils2plot, tok);

for i = idx2plot(:)'

  ref = zeros(tok.nc,1);
  ref(i) = stepsize;
  
  out = run_coilsim(G, K, circ.vmax, circ.vmin, tmax, ref);
  
  out.e = resample(out.e, t);
  out.x = resample(out.x, t);
  out.u = resample(out.u, t);

  e(:,:,i) = out.e.Data;
  x(:,:,i) = out.x.Data;
  u(:,:,i) = out.u.Data;
end

cc_analysis_plots(t, x, u, L, tok, coils2plot);

























