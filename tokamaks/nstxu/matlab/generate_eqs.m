clear all; clc; close all

saveit = 1;
savefn = './eq/eqs3333.mat';

% load circ and tok
circ = nstxu_circ();
tok = load('nstxu_obj_2016_GSgrid33x33_npp4x4.mat').tok_data_struct;
tok = connect_tok(tok, circ);


% this will hold all the equilibria we create
eqs = {};  


% Create 3 lower single null equilibria at different elongation
for elong = [1.7 1.9 2.1]

  [spec, init, config] = gsdesign_nstxu_easyload('LSN', tok, circ);

  r = spec.targets.rsep;
  z = spec.targets.zsep;
  s = shape_analysis(r,z);
  s.elong = elong;              % modify the elongation
  [r,z] = shape_edit(r,z,s);     
  spec.targets.rsep = r;
  spec.targets.zsep = z;

  eqs{end+1} = gsdesign(spec, init, config);
end



% Create 3 lower-single null equilibria at different triangularity
for dtri = [-0.2 -0.1 0.2]

  [spec, init, config] = gsdesign_nstxu_easyload('LSN', tok, circ);

  r = spec.targets.rsep;
  z = spec.targets.zsep;
  s = shape_analysis(r,z);
  s.triu = s.triu + dtri;         % modify the triangularity
  s.tril = s.tril + dtri;         
  [r,z] = shape_edit(r,z,s);     
  spec.targets.rsep = r;
  spec.targets.zsep = z;

  eqs{end+1} = gsdesign(spec, init, config);
end



% Create 3 limited equilibria at different elongation
for elong = [1.4 1.6 1.8]

  [spec, init, config] = gsdesign_nstxu_easyload('LIM', tok, circ);

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










