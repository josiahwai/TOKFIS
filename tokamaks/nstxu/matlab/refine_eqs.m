clear all; clc; close all

% convert 33x33 grid equilibria to 65x65 grid

saveit = 1;
savefn = './eq/eqs6565.mat';
loadfn = './eq/eqs3333.mat';


% load stuff
circ = nstxu_circ();
tok = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
tok = connect_tok(tok, circ);
eqs = load(loadfn).eqs;


for i = 1:length(eqs)

  [spec, config] = gsdesign_recreate_init(eqs{i}, tok);
  
  spec.mxiter = 12;
  spec.weights.sep(:) = 1e3;
  spec.targets.ic = eqs{i}.ic;
  spec.weights.ic = ones(size(spec.targets.ic)) * 1e-2;

  eq = gsdesign(spec, eqs{i}, config);
  eq = rmfield(eq, {'p','b','r','d','source'});  % these take up a lot of space
  eqs{i} = eq;

end

% save stuff
if saveit
  save(savefn, 'eqs')
end










