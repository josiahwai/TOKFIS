
clc; close all

iuse = [4 14 15 25];
model_order = 20;

G = {};
ip = {};
for i = 1:length(iuse)  
  sys = models{iuse(i)};
  ip{i} = sys.eq.cpasma;
  G{i} = ss(sys.amat, sys.bmat, sys.dzcurdx, 0);
  G{i} = reduce(G{i}, model_order);  
end

s = tf('s');

kp = 1e-3;
kd = 8e-5;
ki = 1e-3;
delay = 1e-3;

ssdelay = exp(-delay*s);
ssdelay = pade(ssdelay, 4);

K = struct;
for dum = {'PF1', 'PF2', 'PF3'}
  x = dum{:};

  K.(x) = tunablePID(x, 'pid');
  K.(x).Kp.Value = kp;
  K.(x).Kp.Minimum = -kp*2;
  K.(x).Kp.Maximum = kp*10;
  K.(x).Kp.Free = 1;
  K.(x).Ki.Value = ki * 0;
  K.(x).Ki.Minimum = -ki*2;
  K.(x).Ki.Maximum = ki * 10;
  K.(x).Ki.Free = 0;
  K.(x).Kd.Value = kd;
  K.(x).Kd.Minimum = -kd*2;
  K.(x).Kd.Maximum = kd*10;
  K.(x).Kd.Free = 1;
  K.(x).Tf.Value = 1e-6;
  K.(x).Tf.Free = 0;
end

vec1 = [0 1 0 0 0 0 0 -1]';
vec2 = [0 0 1 0 0 0 -1 0]';
vec3 = [0 0 0 1 0 -1 0 0]';


CL0 = {};
for i = 1:length(G)
  Kvs = ip{i} * ssdelay * (K.PF1*vec1 + K.PF2*vec2 + K.PF3*vec3);

  apnam1 = ['X' num2str(i)];
  ap1 = AnalysisPoint(apnam1);  % analysis point on feedback loop

  apnam2 = ['E' num2str(i)];
  ap2 = AnalysisPoint(apnam2);

  CL0{i} = feedback(G{i}*Kvs*ap2, ap1);
end
CL0 = vertcat(CL0{:});



%%
clc
close all

CL0.InputName = 'r';
CL0.OutputName = {'z1', 'z2', 'z3', 'z4'};

ny = size(CL0,1);
hardreqs = {};
softreqs = {};




% % Tracking goals
% responsetime = 0.15;
% dcerror      = 0.1;
% for i = 1:ny
%   softreqs{end+1} = TuningGoal.Tracking('r', CL0.OutputName{i}, responsetime, dcerror);
% end


% % Overshoot
% maxpercent = 100;
% for i = 1:ny
%   hardreqs{end+1} = TuningGoal.Overshoot('r', CL0.OutputName{i}, maxpercent);
% end


% % Gain
% gainvalue = 2.2;
% for i = 1:ny
%   hardreqs{end+1} = TuningGoal.Gain('r', CL0.OutputName{i}, gainvalue);
% end


% % Margins
% gainmargin = 1;
% phasemargin = 15;
% for i = 1:ny
%   apnam = ['X' num2str(i)];
%   softreqs{end+1} = TuningGoal.Margins(apnam, gainmargin,phasemargin);
% end


% Frequency-based gain
hfgain = 3;
lfgain = 0.05;
wc = 60;      % rad/s
W = makeweight(lfgain, [wc 1], hfgain);

% hold on
% bodemag(W)

for i = 1:ny
  ap = ['E' num2str(i)];
  hardreqs{end+1} = TuningGoal.Gain('r', ap, W);
end


softreqs = [softreqs{:}];
hardreqs = [hardreqs{:}];

opt = systuneOptions;
opt.SoftTol = 1e-4;
opt.RandomStart = 2;
opt.MaxIter = 100;

[CL,fSoft] = systune(CL0, softreqs, hardreqs, opt);



[y0,t0] = step(CL0, 0.5);
[y,t] = step(CL, 0.5);
figure
hold on
% plot(t0, y0, '--r')
plot(t, y)
axis([0 0.5 -0.5 2.5])


K = struct;
for dum = {'PF1', 'PF2', 'PF3'}
  x = dum{:};
  K.(x) = getBlockValue(CL, x);
end


%%

G = {};
Gcl = {};

for i = 1:ny

  sys = models{iuse(i)};
  ip{i} = sys.eq.cpasma;
  G{i} = ss(sys.amat, sys.bmat, eye(size(sys.amat)), 0);

  Kvs = ip{i} * ssdelay * (K.PF1*vec1 + K.PF2*vec2 + K.PF3*vec3); 
  Gcl{i} = feedback(G{i}, Kvs * sys.dzcurdx);

end

%%


% 
% figure
% hold on
% for i = 1:ny
%   T = CL(i);
%   S = 1 - T;
%   L = inv(S) - 1;
% 
%   nyquist(L)
%   axis([-2 0 -2 0.2])
% 
%   plot([-2 0], [-1 0], 'k', 'linewidth', 2)
%    plot([-2 0], [-2 0], 'k', 'linewidth', 2)
% end












































