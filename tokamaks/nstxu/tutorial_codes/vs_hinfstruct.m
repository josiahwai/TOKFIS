% Structured Hinf to tune the RZIp PID controllers
clear all; clc; close all

%% load models
models = load('../mat2/models_gsupdate.mat').models;

% iuse = [100 200]';
% iuse = [6, 25, 30, 58, 59, 65, 77, 89, 125, 150, 245, 255, 275]';
load('iuse');
N = length(iuse);

G = {};
for k = 1:N
  
  ii = iuse(k);
  ip = models.ip(ii);

  % define model
  A = squeeze(models.A(ii,:,:));
  B = squeeze(models.B(ii,:,:));
  Crzip = squeeze(models.Crzip(ii,:,:));
  Crzip = diag([ip ip 1]) * Crzip;        % r and z are scaled propto ip
  nx = size(A,1);
  nu = size(B,2);

  P = ss(A, B, Crzip, 0);
  
  % normalize
  Du = 1;                              
  Dx = diag( [ones(1,nx-1) 1e3]); 
  Dy = diag([1 1 1e3]);              
  P = scalesys(P,Du,Dx,Dy);
  P = P(2,:);

  G{k} = P;  % transfer fun from voltage to zcur
 
  
%   [bsys,g,T,Ti] = balreal(P);
%   n = 25;
%   elim = false(1,length(g));
%   elim(n+1:end) = 1;
%   P = modred(bsys,elim);
%   G{k} = P;
  
end

%% Tunable parameters

kpf1 = tunablePID('kpf1','pd'); 
kpf1.Kp.Value = 1e-5;
kpf1.Kp.Minimum = 0;
kpf1.Kp.Maximum = 1e-2;
kpf1.Kd.Value = 1e-5;
kpf1.Kd.Minimum = 0;
kpf1.Kd.Maximum = 1e-2;
kpf1.Tf.Value = 1e-6;
kpf1.Tf.Free = 0;

kpf2 = tunablePID('kpf2','pd'); 
kpf2.Kp.Value = 1e-5;
kpf2.Kp.Minimum = 0;
kpf2.Kp.Maximum = 1e-2;
kpf2.Kd.Value = 1e-5;
kpf2.Kd.Minimum = 0;
kpf2.Kd.Maximum = 1e-2;
kpf2.Tf.Value = 1e-6;
kpf2.Tf.Free = 0;

kpf3 = tunablePID('kpf3','pd'); 
kpf3.Kp.Value = 2e-4;
kpf3.Kp.Minimum = 0;
kpf3.Kp.Maximum = 1e-2;
kpf3.Kd.Value = 1e-4;
kpf3.Kd.Minimum = 0;
kpf3.Kd.Maximum = 1e-2;
kpf3.Tf.Value = 1e-6;
kpf3.Tf.Free = 0;

Kz0 = kpf1 * [0 1 0 0 0 0 0 -1]' + ...
      kpf2 * [0 0 1 0 0 0 -1 0]' + ...
      kpf3 * [0 0 0 1 0 -1 0 0]';

% Hinf performance weights
W1 = makeweight(60, 40, 0.1);
W2 = [];
W3 = 0.8;
% bodemag(W1); hold on

CL0 = [];
for i = 1:N
  genP = augw(G{i}, W1, W2, W3);
  CL0i = lft(genP, Kz0);
  CL0 = blkdiag(CL0, CL0i);
end




%% structured Hinf synthesis
rng('default')
opt = hinfstructOptions('Display','iter','RandomStart',0,'MaxIter',20);
% CL = hinfstruct(genP, zvec*Kz0,opt);
CL = hinfstruct(CL0, opt);

% load tuned parameters
kpf1 = getBlockValue(CL, 'kpf1');
kpf2 = getBlockValue(CL, 'kpf2');
kpf3 = getBlockValue(CL, 'kpf3');

Kz = kpf1 * [0 1 0 0 0 0 0 -1]' + ...
      kpf2 * [0 0 1 0 0 0 -1 0]' + ...
      kpf3 * [0 0 0 1 0 -1 0 0]';


%% plot results

% s = tf('s');
% Kz0 = 4e-3 + 2e-4*s;
% Kz = Kz0;
% zvec = zvec0;

T = {}; L = {}; T0 = {};  L0 = {}; re = []; im = []; y = [];  y0 = [];

for i = 1:N
  L{i} = G{i}*Kz;
  T{i} = feedback(L{i}, 1);
  
  L0{i} = ss(G{i}*Kz0);
  T0{i} = feedback(L0{i}, 1);
  
  % step responses
  t = linspace(0, 0.2, 100);
  y(i,:) = step(T{i}, t);
  y0(i,:) = step(T0{i}, t);
  
  % nyquists
  w = logspace(3, 0, 30);
  w = [flip(-w) w];
  [re(i,:), im(i,:)] = nyquist(L{i}, w);  
  [re0(i,:), im0(i,:)] = nyquist(L0{i}, w);  
end
  
f = figure;
f.Position = [402 323 702 290];
subplot(121)
hold on
grid on
plot(t,y0,'color', [1 1 1] * 0.7) 
subplot(122)
hold on
grid on
plot(t,y,'color', [1 1 1] * 0.7) 

f = figure;
f.Position = [489 397 755 278];
subplot(121)
hold on
grid on
plot(re0', im0', 'color', [1 1 1] * 0.7)
scatter(-1, 0, 20, 'k', 'filled')
axis([-3 1 -2 2])
subplot(122)
hold on
grid on
plot(re', im', 'color', [1 1 1] * 0.7)
scatter(-1, 0, 20, 'k', 'filled')
axis([-3 1 -2 2])





















































