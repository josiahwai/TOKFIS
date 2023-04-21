clear all; clc; close all

% convert 33x33 grid equilibria to 65x65 or 129x129 grid

saveit = 1;
savefn = './eq/eqs129.mat';
loadfn = './eq/eqs3333.mat';


% load stuff
circ = sparc_circ();
tok = load('sparc_obj_129.mat').tok_data_struct;
tok = connect_tok(tok, circ);
eqs = load(loadfn).eqs;


for i = 1:length(eqs)
  i
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










