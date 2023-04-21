clear all; clc; close all
load('x')

gbr = [tok.gbr2c tok.gbr2v];
gbz = [tok.gbz2c tok.gbz2v];
mpx = [tok.mpc tok.mpv];

% vacuum fields
x = [eq.ic; eq.iv];
br = reshape(gbr * x, tok.nz, tok.nr);
bz = reshape(gbz * x, tok.nz, tok.nr);

in = inpolygon(tok.rgg, tok.zgg, tok.limdata(2,:), tok.limdata(1,:));
br(~in) = nan;
bz(~in) = nan;

dr = mean(diff(tok.rg));
dz = mean(diff(tok.zg));
[dbrdr, dbrdz] = gradient(br, dr, dz);
[dbzdr, dbzdz] = gradient(bz, dr, dz);


dbrdz = dbrdz * 1e2;  % [T/m] --> [Gs/cm]
dbzdr = dbzdr * 1e2;  
dbrdr = dbrdr * 1e2;
dbzdz = dbzdz * 1e2;

figure
contourf(tok.rg, tok.zg, dbrdz, 50);
colorbar; plot_lim(tok);

i = tok.rgg < 2.2 & tok.rgg > 1.4 & tok.zgg < 0.4 & tok.zgg > -0.4;
mean(abs(dbrdz(i)))


figure
contourf(tok.rg, tok.zg, dbzdr, 50);
colorbar; plot_lim(tok);

i = tok.rgg < 2.2 & tok.rgg > 1.4 & tok.zgg < 0.4 & tok.zgg > -0.4;
mean(abs(dbzdr(i)))




% br(~in) = nan;
% figure
% contourf(tok.rg, tok.zg, br)
% colorbar



% gbr = [tok.gbr2c tok.gbr2v];
% gbz = [tok.gbz2c tok.gbz2v];
% mpx = [tok.mpc tok.mpv];
% 
% % vacuum fields
% x = [eq.ic; eq.iv];
% brvac = reshape(gbr * x, nz, nr);
% bzvac = reshape(gbz * x, nz, nr);
% 
% % current distribution and gradient
% J = eq.pcurrt;
% ip = sum(J(:));
% rc = sum(J(:).*rgg(:)) / ip;
% zc = sum(J(:).*zgg(:)) / ip;
% [djdr, djdz] = gradient(J, dr, dz);











