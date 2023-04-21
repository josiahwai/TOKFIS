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

% compute plasma coupling
dpsipladx = [tok.mpc tok.mpv] * vec(1:end-1);
dpsipladx = reshape(dpsipladx, tok.nz, tok.nr);


% plots
figure

subplot(2,2,1)
labels = categorical(tok.ccnames);
bar(labels, vec(1:tok.nc))

subplot(2,2,3)
bar(vec)
xlabel('Vessel element')

subplot(2,2,[2 4])
title(['Gamma = ' num2str(val)])
hold on
contourf(tok.rg, tok.zg, dpsipladx, 20)
colorbar 
plot(eq.rbbbs, eq.zbbbs, 'w', 'linewidth', 1.5)
plot_lim(tok)













