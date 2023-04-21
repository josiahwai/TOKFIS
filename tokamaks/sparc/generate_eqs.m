clear all; clc; close all

saveit = 0;
savefn = './eq/eqs3333.mat';

% load circ and tok
circ = sparc_circ();
tok = load('sparc_obj_3333.mat').tok_data_struct;
tok = connect_tok(tok, circ);


% this will hold all the equilibria we create
eqs = {};  

 
% Create 3 lower single null equilibria at different elongation
for elong = [1.8 1.9 2]

  [spec, init, config] = gsdesign_sparc_easyload('LSN', tok, circ);

  r = spec.targets.rsep;
  z = spec.targets.zsep;
  s = shape_analysis(r,z);
  s.elong = elong;              % modify the elongation
  [r,z] = shape_edit(r,z,s);     
  spec.targets.rsep = r;
  spec.targets.zsep = z;

  eqs{end+1} = gsdesign(spec, init, config);
end



% Create 3 equilibria with different li
for li = [0.7 1 1.3]

  [spec, init, config] = gsdesign_sparc_easyload('LSN', tok, circ);

  spec.targets.li = li;
  spec.weights.li = 1e3;

  eqs{end+1} = gsdesign(spec, init, config);
end



% Create 3 limited equilibria at different elongation
for elong = [1.4 1.6 1.8]

  [spec, init, config] = gsdesign_sparc_easyload('LIM', tok, circ);

  r = spec.targets.rsep;
  z = spec.targets.zsep;
  s = shape_analysis(r,z);
  s.elong = elong;              % modify the elongation
  [r,z] = shape_edit(r,z,s);     
  spec.targets.rsep = r;
  spec.targets.zsep = z;

  eqs{end+1} = gsdesign(spec, init, config);
end



% save eqs
% remove stuff that takes up a lot of space
if saveit
  for i = 1:length(eqs)
    fds2remove = {'p','b','r','d','source'};
    try
      eqs{i} = rmfield(eqs{i}, fds2remove);
    catch
    end
  end
  
  save(savefn, 'eqs')
end










