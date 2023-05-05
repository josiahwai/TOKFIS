function Actuators = nstxu_Actuators()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   Actuators = nstxu_Actuators()
%
%  PURPOSE: Setup for block Actuators in nstxu.slx
%
%  INPUTS: s, simulation settings
%          Tokamak, output from nstxu_Tokamak.m
%
%  OUTPUT: Actuators, data for subsystem Actuators in nstxu.slx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%
%  WRITTEN BY:  Anders Welander ON 2021-01-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PS.names = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'}';
nPS = length(PS.names);

PS.curdir = [0 1 1 0 -1 0 1 1];   % 0=bipolar, +/-1 = positive/negative current only 
PS.vmax = [4048 1012 2024 2024 3036 2024 2024 1012]';
PS.vmin = -PS.vmax;

frequency = 65;
PS.slewratemax = 2*frequency*PS.vmax;
PS.slewratemin = 2*frequency*PS.vmin;
PS.latency = 100e-6;

Actuators.PS = PS;
Actuators.nPS = nPS;





