clear all; clc; close all

% eqs_fn = './eq/eqs33_test.mat';
% eqs_fn = './eq/eqs3333.mat';
eqs_fn = './eq/eqs6565.mat';
% eqs_fn = './eq/eqs129.mat';

% tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';
% tok_fn = 'nstxu_obj_config2016_6565.mat';
% tok_fn = 'kstar_obj_mks_struct_2015_3333.mat';
% tok_fn = 'kstar_obj_mks_struct_2015_6565.mat';
% tok_fn = 'kstar_obj_mks_struct_2017_3333.mat';
% tok_fn = 'sparc_obj_3333.mat';
tok_fn = 'sparc_obj_6565.mat';
% tok_fn = 'sparc_obj_129.mat';

% circ = nstxu_circ();
% circ = kstar_circ();
circ = sparc_circ();

% tok = tok_data_struct;
% A = -inv(tok.mvv)*diag(tok.resv); e = esort(eig(A)); e(1:5)

% load stuff
tok = load(tok_fn).tok_data_struct;
eqs = load(eqs_fn).eqs;

% tok = rmfield(tok, {'nrl','rldata','rlsignals','rlnames'});

% connect coils-->circuits
tok = connect_tok(tok, circ);          

%%
eq = eqs{2};

% plasma_model = {'rzrig', 'rig', 'gspert'};
plasma_model = 'all';
iplcirc = 1; 
Rp = 5e-7;
models = response_models(eq, tok, iplcirc, Rp, plasma_model);

e = [];
e(:,1) = esort(eig(models.rzrig.amat));    
e(:,2) = esort(eig(models.gspert.amat));   
e(:,3) = esort(eig(models.gsupdate.amat)); 
e(:,4) = esort(eig(models.rig.amat));      
e(:,5) = esort(eig(models.vacuum.amat));   

real(e(1:5,:))


%%

if 0
  close all
  k = 1;
  
  visualize_eig_vec(models.rzrig.amat, eq, tok, k)
  sgtitle('rzrig')
  
  visualize_eig_vec(models.rig.amat, eq, tok, k)
  sgtitle('rig')
  
  visualize_eig_vec(models.gspert.amat, eq, tok, k)
  sgtitle('gspert')
  
  % visualize_eig_vec(models.gsupdate.amat, eq, tok, k)
  % sgtitle('gsupdate')
  
  % visualize_eig_vec(models.vacuum.amat, eq, tok, k)
  % sgtitle('vacuum')

end




















































