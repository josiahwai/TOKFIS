clear all; clc; close all

% tok_fn = 'nstxu_obj_2016_GSgrid33x33_npp4x4.mat';

tok_fn = 'nstxu_obj_config2016_6565.mat';
eqs_fn = './eq/eqs6565.mat';

% LOAD SOME EQUILIBRIA AND GEOMETRY
% (to create different equilibria, see generate_eqs.m)
tok = load(tok_fn).tok_data_struct;
eqs = load(eqs_fn).eqs;
eq = eqs{3};

% CONNECT COILS-->CIRCUITS
circ = nstxu_circ();
tok = connect_tok(tok, circ);          


% BUILD THE VACUUM MODEL
iplcirc = 1;
Rp = 5e-7;
plasma_model = 'rzrig';
% sys = response_models(eq, tok, iplcirc, Rp, plasma_model)



mpc = tok.mpc;
mps = [tok.mpc tok.mpv];
mss = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
dz = mean(diff(tok.zg));
dr = mean(diff(tok.rg));
mpc = reshape(mpc, tok.nz, tok.nr, []);
mps = reshape(mps, tok.nz, tok.nr, []);
[~, dmpcdz] = gradient(mpc, dr, dz, 1);
[~, dmpsdz] = gradient(mps, dr, dz, 1);
[~, d2mpcdz2] = gradient(dmpcdz, dr, dz, 1);

dmpcdz = reshape(dmpcdz, tok.nz*tok.nr, []);
dmpsdz = reshape(dmpsdz, tok.nz*tok.nr, []);
d2mpcdz2 = reshape(d2mpcdz2, tok.nz*tok.nr, []);


mpp = unwrap_mpp(tok.mpp, tok.nz, tok.nr);
mpp = reshape(mpp, tok.nz, tok.nr, []);
[~, dmppdz] = gradient(mpp, dr, dz, 1);
[~, d2mppdz2] = gradient(dmppdz, dr, dz, 1);
d2mppdz2 = reshape(d2mppdz2, tok.nz*tok.nr, []);

%%
ic = eq.ic;
% ic([4 6]) = ic([4 6]) + * 1000;
ic(4) = ic(4) + 100;

f = (eq.pcurrt(:)' * dmpsdz * inv(mss) * dmpsdz' * eq.pcurrt(:)) / ...
  (eq.pcurrt(:)' * d2mpcdz2 * ic)


%%
n = 100;

c1s = linspace(0,1,n);
c3s = linspace(0,1,n);

num = eq.pcurrt(:)' * dmpsdz * inv(mss) * dmpsdz' * eq.pcurrt(:);
f0 = num/ (eq.pcurrt(:)' * d2mpcdz2 * eq.ic);
% for i = 1:n
%   for j = 1:n
%     c1 = c1s(i);
%     c3 = c3s(j);
%     vec = normc([c1; c3]) * 100;
% 
%     ic = eq.ic;
%     ic([2 4]) = ic([2 4]) + vec;
% 
%     f(i,j) = num/ (eq.pcurrt(:)' * d2mpcdz2 * ic);
%   end
% end

f = [];
for i = 1:n
  c1 = c1s(i);
  c3 = sqrt(1-c1^2);
  vec = -[c1; c3] * 100;

  ic = eq.ic;
  ic([2 4]) = ic([2 4]) + vec;
  % ic([6 8]) = ic([6 8]) + vec;

  f(i) = num/ (eq.pcurrt(:)' * d2mpcdz2 * ic);
end


hold on
yline(f0, 'k')
plot(c1s, f)





%%
num = eq.pcurrt(:)' * dmpsdz * inv(mss) * dmpsdz' * eq.pcurrt(:);
    
f = [];
for i = 1:tok.nz
  for j = 1:tok.nr
    
    sz = [tok.nz tok.nr];
    k = sub2ind(sz, i, j);
    x = -d2mppdz2(:,k) * 1000;

    
    
    den = eq.pcurrt(:)' * (d2mpcdz2 * ic + x);

    f(i,j) = num/den;

  end
end

clf
plot_lim(tok)
contour(tok.rg, tok.zg, f, 20)
plot(eq.rbbbs, eq.zbbbs, 'k')
colorbar



%%
n = 20;
th = linspace(0,pi/2,n);
phi = linspace(0,pi/2,n);
[th,phi] = meshgrid(th,phi);
r = ones(size(th)) * 1;
[x,y,z] = sph2cart(th, phi, r);
x = x(:);
y = y(:);
z = z(:);


f = [];
for i = 1:length(x)


  vec = [x(i) y(i) z(i)]';
  % vec = [x(i) z(i) y(i)]';
  % vec = [z(i) y(i) x(i)]';

  ic = eq.ic;
  ic([2 3 4]) = ic([2 3 4]) - vec*100;

   f(i) = num/ (eq.pcurrt(:)' * d2mpcdz2 * ic);

end


%Create regular grid across data space
n = 100;
[X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));

%create contour plot
F = griddata(x,y,f,X,Y);
contourf(X,Y,F,20)
colorbar;
xlabel('c1')
ylabel('c2')

%%



c1s = linspace(0,1,n);
c3s = linspace(0,1,n);

num = eq.pcurrt(:)' * dmpsdz * inv(mss) * dmpsdz' * eq.pcurrt(:);
f0 = num/ (eq.pcurrt(:)' * d2mpcdz2 * eq.ic);
% for i = 1:n
%   for j = 1:n
%     c1 = c1s(i);
%     c3 = c3s(j);
%     vec = normc([c1; c3]) * 100;
% 
%     ic = eq.ic;
%     ic([2 4]) = ic([2 4]) + vec;
% 
%     f(i,j) = num/ (eq.pcurrt(:)' * d2mpcdz2 * ic);
%   end
% end

f = [];
for i = 1:n
  c1 = c1s(i);
  c3 = sqrt(1-c1^2);
  vec = -[c1; c3] * 100;

  ic = eq.ic;
  ic([2 4]) = ic([2 4]) + vec;
  % ic([6 8]) = ic([6 8]) + vec;

  f(i) = num/ (eq.pcurrt(:)' * d2mpcdz2 * ic);
end


hold on
yline(f0, 'k')
plot(c1s, f)


