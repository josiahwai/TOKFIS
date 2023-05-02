% This example shows how to use systune to design a PID controller 
%
% The controller searches for kp, ki, kd coefficients for PF1, PF2, PF3
% coil sets to meet the target design. 

function [Kz, K] = vs_optimization(models)

% define models
% transfer function from voltage to zcur
model_order = 20;
G = {};
ip = {};
for i = 1:length(models)
  sys = models{i};
  ip{i} = sys.eq.cpasma;
  G{i} = ss(sys.amat, sys.bmat, sys.dzcurdx, 0);
  G{i} = reduce(G{i}, model_order);  
end



% initial values for controller
kp = 1e-3;       
kd = 8e-5;
ki = 1e-3;
delay = 1e-3;    % controller delay [s]

% approximate time delay using pade
s = tf('s');
ssdelay = exp(-delay*s);
ssdelay = pade(ssdelay, 4);


% initialize controller by defining tunable PID blocks
K = struct;
for dum = {'PF1', 'PF2', 'PF3'}
  x = dum{:};

  K.(x) = tunablePID(x, 'pid');
  K.(x).Kp.Value = kp;
  K.(x).Kp.Minimum = -kp*2;
  K.(x).Kp.Maximum = kp*10;
  K.(x).Kp.Free = 1;
  K.(x).Ki.Value = ki;
  K.(x).Ki.Minimum = -ki*2;
  K.(x).Ki.Maximum = ki * 10;
  K.(x).Ki.Free = 1;
  K.(x).Kd.Value = kd;
  K.(x).Kd.Minimum = -kd*2;
  K.(x).Kd.Maximum = kd*10;
  K.(x).Kd.Free = 1;
  K.(x).Tf.Value = 1e-6;
  K.(x).Tf.Free = 0;
end


vec1 = [0 1 0 0 0 0 0 -1]';  % vector representation of (+PF1U, -PF1L)
vec2 = [0 0 1 0 0 0 -1 0]';  % vector representation of (+PF2U, -PF2L)
vec3 = [0 0 0 1 0 -1 0 0]';  % vector representation of (+PF3U, -PF3L)


% close the vertical control loop 
CL0 = {};
for i = 1:length(G)

  Kvs = ip{i} * ssdelay * (K.PF1*vec1 + K.PF2*vec2 + K.PF3*vec3);   % representation of vertical controller

  apnam1 = ['X' num2str(i)];    % analysis point on feedback loop, at output y location
  ap1 = AnalysisPoint(apnam1);  

  apnam2 = ['E' num2str(i)];    % analysis point at error location
  ap2 = AnalysisPoint(apnam2);

  CL0{i} = feedback(G{i}*Kvs*ap2, ap1);   % close the loop
end

CL0 = vertcat(CL0{:});  % concatenate all models together



%% Define tuning and tracking goals

ny = size(CL0,1);       % number of outputs
CL0.InputName = 'r';    % define input names

for i = 1:ny
  CL0.OutputName{i} = ['z' num2str(i)];   % define outputnames
end

% these will hold all the hard and soft requirements for systune
hardreqs = {};
softreqs = {};


% Note: these tuning goals are just for reference, may not actually specify
% a good design. Modify the if statements to add or remove goals. Also,
% user may want to change whether a goal belongs to softreqs or hardreqs.
% In my experience the 'simpler' goals like overshoot don't always create a good
% response and the frequency-based weighting functions do beter. 

% Tracking goals
if 0
  responsetime = 0.15;
  dcerror      = 0.1;
  for i = 1:ny
    softreqs{end+1} = TuningGoal.Tracking('r', CL0.OutputName{i}, responsetime, dcerror);
  end
end


% Overshoot goal
if 0
  maxpercent = 100;
  for i = 1:ny
    hardreqs{end+1} = TuningGoal.Overshoot('r', CL0.OutputName{i}, maxpercent);
  end
end


% Gain goals
if 0
  gainvalue = 2.2;
  for i = 1:ny
    hardreqs{end+1} = TuningGoal.Gain('r', CL0.OutputName{i}, gainvalue);
  end
end

% Margins
if 0
  gainmargin = 1;
  phasemargin = 15;
  for i = 1:ny
    apnam = ['X' num2str(i)];
    softreqs{end+1} = TuningGoal.Margins(apnam, gainmargin,phasemargin);
  end
end


% Frequency-based gain on the sensitivity function (reference-->error)
if 1
  hfgain = 3;
  lfgain = 0.05;
  wc = 60;      % rad/s
  W = makeweight(lfgain, [wc 1], hfgain);
  
  for i = 1:ny
    ap = ['E' num2str(i)];
    hardreqs{end+1} = TuningGoal.Gain('r', ap, W);
  end
end


% Combine all soft and hard requirements and synthesize the controller
softreqs = [softreqs{:}];
hardreqs = [hardreqs{:}];

opt = systuneOptions;
opt.SoftTol = 1e-4;
opt.RandomStart = 2;
opt.MaxIter = 100;

[CL,fSoft] = systune(CL0, softreqs, hardreqs, opt);


% Plot step response of system
[y,t] = step(CL, 0.5);
figure
hold on
plot(t, y)
axis([0 0.5 -0.5 2.5])
xlabel('Time [s]')
ylabel('Z position')
title('Step response of synthesized controller', 'fontsize', 14)


% Read controller values
K = struct;
for dum = {'PF1', 'PF2', 'PF3'}
  x = dum{:};
  K.(x) = getBlockValue(CL, x);
end


Kz = (K.PF1*vec1 + K.PF2*vec2 + K.PF3*vec3);










































