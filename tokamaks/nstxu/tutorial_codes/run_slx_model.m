% Helper function for programmatically calling and running a simulink slx
% model with fast restart.

function out = run_slx_model(filename, varargin)

% extract filename without extension
[pth, nam, ext] = fileparts(filename);


% create Simulation Input object
simIn = Simulink.SimulationInput(nam);

% add all variables to simulink 
for i = 1:length(varargin)
  fd = inputname(i+1);
  val = varargin{i};
  simIn = setVariable(simIn, fd, val);
end

% run simulation with fast restart 
out = sim(simIn, 'UseFastRestart', 'on');












