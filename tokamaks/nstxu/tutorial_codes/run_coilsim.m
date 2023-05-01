function out = run_coilsim(G, K, vmax, vmin, tmax, ref)

simIn = Simulink.SimulationInput('coilsim_slx');

G = ss(G);
K = ss(K);
Cc = eye(size(K,2), size(G,1));
s = variables2struct(G,K,Cc,vmax,vmin,tmax,ref);

for x = fields(s)'
  fd = x{:};
  simIn = setVariable(simIn, fd, s.(fd));
end
 
out = sim(simIn, 'UseFastRestart', 'on');










