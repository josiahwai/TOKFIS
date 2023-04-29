% A   - dynamics matrix
% tok - tok_data_struct geometry
% k   - index of eigenvector to visualise (1 = most unstable) 
%
function visualize_eig_vec(A, eq, tok, k)


% compute eigenvalues
[V,D] = eig(A);

% sort eigenvalues
[d,i] = esort(diag(D));
V = V(:,i);

% index a specific eigenvalue
val = d(k);
vec = V(:,k);
ic = vec(1:tok.nc);
iv = vec(tok.nc+1:tok.nc+tok.nv);

% compute plasma coupling
dpsipladx = tok.mpc*ic + tok.mpv*iv;
dpsipladx = reshape(dpsipladx, tok.nz, tok.nr);


% plots
figure

subplot(2,2,1)
labels = categorical(tok.ccnames);
bar(labels, ic)

subplot(2,2,3)
bar(iv)
xlabel('Vessel element')

subplot(2,2,[2 4])
title(['Gamma = ' num2str(val)])
hold on
contourf(tok.rg, tok.zg, dpsipladx, 20)
colorbar 
plot(eq.rbbbs, eq.zbbbs, 'w', 'linewidth', 1.5)
plot_lim(tok)













